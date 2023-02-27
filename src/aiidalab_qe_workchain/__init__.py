# AiiDA imports.
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, if_
from aiida.orm import CalcJobNode, QueryBuilder, WorkChainNode
from aiida.plugins import DataFactory, WorkflowFactory

# AiiDA Quantum ESPRESSO plugin inputs.
from aiida_quantumespresso.common.types import RelaxType
from aiida_quantumespresso.utils.mapping import prepare_process_inputs

# Data objects and work chains.
PwRelaxWorkChain = WorkflowFactory("quantumespresso.pw.relax")
PwBandsWorkChain = WorkflowFactory("quantumespresso.pw.bands")
PdosWorkChain = WorkflowFactory("quantumespresso.pdos")
XspectraCoreWorkChain = WorkflowFactory("quantumespresso.xspectra.core")
XpsWorkChain = WorkflowFactory("quantumespresso.xps")

Bool = DataFactory("core.bool")
Float = DataFactory("core.float")
Dict = DataFactory("core.dict")
Str = DataFactory("core.str")
XyData = DataFactory("core.array.xy")
StructureData = DataFactory("core.structure")
BandsData = DataFactory("core.array.bands")
Orbital = DataFactory("core.orbital")


def get_pseudo(group, element):
    pseudos = {pseudo.element: pseudo for pseudo in group.nodes}
    return pseudos.get(element, None)


class QeAppWorkChain(WorkChain):
    """WorkChain designed to calculate the requested properties in the AiiDAlab Quantum ESPRESSO app."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('structure', valid_type=StructureData,
                   help='The inputs structure.')
        spec.input('clean_workdir', valid_type=Bool, default=lambda: Bool(False),
                   help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.expose_inputs(PwRelaxWorkChain, namespace='relax', exclude=('clean_workdir', 'structure'),
                           namespace_options={'required': False, 'populate_defaults': False,
                                              'help': 'Inputs for the `PwRelaxWorkChain`, if not specified at all, the relaxation step is skipped.'})
        spec.expose_inputs(PwBandsWorkChain, namespace='bands',
                           exclude=('clean_workdir', 'structure', 'relax'),
                           namespace_options={'required': False, 'populate_defaults': False,
                                              'help': 'Inputs for the `PwBandsWorkChain`.'})
        spec.expose_inputs(PdosWorkChain, namespace='pdos',
                           exclude=('clean_workdir', 'structure'),
                           namespace_options={'required': False, 'populate_defaults': False,
                                              'help': 'Inputs for the `PdosWorkChain`.'})
        spec.expose_inputs(XspectraCoreWorkChain, namespace='xspectra',
                           exclude=('clean_workdir', 'structure'),
                           namespace_options={'required': False, 'populate_defaults': False,
                                              'help': 'Inputs for the `XspectraCoreWorkChain`.'})
        spec.expose_inputs(XpsWorkChain, namespace='xps',
                           exclude=('clean_workdir', 'structure', 'relax'),
                           namespace_options={'required': False, 'populate_defaults': False,
                                              'help': 'Inputs for the `XpsWorkChain`.'})
        spec.input(
            'kpoints_distance_override', valid_type=Float, required=False,
            help='Override for the kpoints distance value of all `PwBaseWorkChains` except for the `nscf` calculations.'
        )
        spec.input(
            'degauss_override', valid_type=Float, required=False,
            help='Override for the `degauss` value of all `PwBaseWorkChains` except for the `nscf` calculations.'
        )
        spec.input(
            'smearing_override', valid_type=Str, required=False,
            help='Override for the `smearing` value of all `PwBaseWorkChains` save for the `nscf` calculations.'
        )
        spec.outline(
            cls.setup,
            if_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax
            ),
            if_(cls.should_run_bands)(
                cls.run_bands,
                cls.inspect_bands
            ),
            if_(cls.should_run_pdos)(
                cls.run_pdos,
                cls.inspect_pdos
            ),
            if_(cls.should_run_xspectra)(
                cls.run_xspectra,
                cls.inspect_xspectra
            ),
            if_(cls.should_run_xps)(
                cls.run_xps,
                cls.inspect_xps
            ),
            cls.results
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
                       message='The PwRelaxWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_BANDS',
                       message='The PwBandsWorkChain sub process failed')
        spec.exit_code(404, 'ERROR_SUB_PROCESS_FAILED_PDOS',
                       message='The PdosWorkChain sub process failed')
        spec.exit_code(405, 'ERROR_SUB_PROCESS_FAILED_XPS',
                       message='The XpsWorkChain sub process failed')
        spec.exit_code(406, 'ERROR_SUB_PROCESS_FAILED_XSPECTRA',
                       message='The XspectraCoreWorkChain sub process failed')
        spec.output('structure', valid_type=StructureData, required=False)
        spec.output('band_parameters', valid_type=Dict, required=False)
        spec.output('band_structure', valid_type=BandsData, required=False)
        spec.output('nscf_parameters', valid_type=Dict, required=False)
        spec.output('dos', valid_type=XyData, required=False)
        spec.output('projections', valid_type=Orbital, required=False)
        spec.output('projections_up', valid_type=Orbital, required=False)
        spec.output('projections_down', valid_type=Orbital, required=False)
        spec.output('symmetry_analysis_data', valid_type=Dict, required=False)
        # for xps
        spec.output_namespace('chemical_shifts', valid_type=Dict, dynamic=True, required=False)
        spec.output_namespace('binding_energies', valid_type=Dict, dynamic=True, required=False)
        spec.output_namespace('xps_spectra_cls', valid_type=XyData, dynamic=True, required=False)
        spec.output_namespace('xps_spectra_be', valid_type=XyData, dynamic=True, required=False)
        # yapf: enable

    @classmethod
    def get_builder_from_protocol(
        cls,
        structure,
        pw_code,
        dos_code=None,
        projwfc_code=None,
        xspectra_code=None,
        protocol=None,
        overrides=None,
        relax_type=RelaxType.NONE,
        pseudo_family=None,
        pseudos=None,
        **kwargs,
    ):
        from aiida.orm import Group

        """Return a builder prepopulated with inputs selected according to the chosen protocol."""
        overrides = overrides.get_dict() if overrides else {}
        builder = cls.get_builder()
        builder.structure = structure

        relax_overrides = overrides.get("relax", {})
        if pseudo_family is not None:
            relax_overrides.setdefault("base", {})["pseudo_family"] = pseudo_family

        relax = PwRelaxWorkChain.get_builder_from_protocol(
            code=pw_code,
            structure=structure,
            protocol=protocol,
            overrides=relax_overrides,
            relax_type=relax_type,
            **kwargs,
        )
        relax.pop("structure", None)
        relax.pop("clean_workdir", None)
        relax.pop("base_final_scf", None)
        builder.relax = relax

        bands_overrides = overrides.get("bands", {})
        if pseudo_family is not None:
            bands_overrides.setdefault("scf", {})["pseudo_family"] = pseudo_family
            bands_overrides.setdefault("bands", {})["pseudo_family"] = pseudo_family
        bands = PwBandsWorkChain.get_builder_from_protocol(
            code=pw_code,
            structure=structure,
            protocol=protocol,
            overrides=bands_overrides,
            **kwargs,
        )
        bands.pop("relax")
        bands.pop("structure", None)
        bands.pop("clean_workdir", None)
        builder.bands = bands

        if dos_code is not None and projwfc_code is not None:
            pdos_overrides = overrides.get("pdos", {})
            if pseudo_family is not None:
                pdos_overrides.setdefault("scf", {})["pseudo_family"] = pseudo_family
                pdos_overrides.setdefault("nscf", {})["pseudo_family"] = pseudo_family
            pdos = PdosWorkChain.get_builder_from_protocol(
                pw_code=pw_code,
                dos_code=dos_code,
                projwfc_code=projwfc_code,
                structure=structure,
                protocol=protocol,
                overrides=pdos_overrides,
                **kwargs,
            )
            pdos.pop("structure", None)
            pdos.pop("clean_workdir", None)
            builder.pdos = pdos
        # xas
        if xspectra_code is not None:
            # xspectra_overrides = overrides.get("xspectra", {})
            from aiida import orm

            pseudos = {
                "Si": orm.load_node(323),
                "X": orm.load_node(325),
            }
            core_wfc_data = orm.load_node(324)
            structure = orm.load_node(1352)
            # if pseudo_family is not None:
            #     xspectra_overrides.setdefault("scf", {})[
            #         "pseudo_family"
            #     ] = pseudo_family
            #     xspectra_overrides.setdefault("nscf", {})[
            #         "pseudo_family"
            #     ] = pseudo_family
            xspectra = XspectraCoreWorkChain.get_builder_from_protocol(
                pw_code=pw_code,
                xs_code=xspectra_code,
                structure=structure,
                protocol=protocol,
                core_hole_pseudos=pseudos,
                core_wfc_data=core_wfc_data,
                # overrides=xspectra_overrides,
                # **kwargs,
            )
            # xspectra.pop("structure", None)
            print(xspectra)
            xspectra["xs_plot"]["xspectra"] = xspectra["xs_prod"]["xspectra"]
            xspectra.pop("clean_workdir", None)
            builder.xspectra = xspectra
        # xps
        pseudo_set = Group
        xps_overrides = overrides.get("xps", {})
        # set pseudo for normal elements
        if pseudo_family is not None:
            xps_overrides.setdefault("ch_scf", {})["pseudo_family"] = pseudo_family
        # load pseudo for excited-state and group-state.
        es_pseudo_family = xps_overrides.pop("es_pseudo", "core_hole")
        gs_pseudo_family = xps_overrides.pop("gs_pseudo", "gipaw")
        es_pseudo_family = (
            QueryBuilder()
            .append(pseudo_set, filters={"label": es_pseudo_family})
            .one()[0]
        )
        gs_pseudo_family = (
            QueryBuilder()
            .append(pseudo_set, filters={"label": gs_pseudo_family})
            .one()[0]
        )
        # set pseudo and core hole treatment for element
        pseudos = {}
        core_hole_treatments = {}
        core_hole_treatment = xps_overrides.pop("core_hole_treatment", "full")
        elements_list = xps_overrides.pop("elements_list", None)
        if not elements_list:
            elements_list = [kind.symbol for kind in structure.kinds]
        for element in elements_list:
            es_pseudo = get_pseudo(es_pseudo_family, element)
            gs_pseudo = get_pseudo(gs_pseudo_family, element)
            if es_pseudo is not None and gs_pseudo is not None:
                pseudos[element] = {
                    "core_hole": es_pseudo,
                    "gipaw": gs_pseudo,
                }
            core_hole_treatments[element] = core_hole_treatment
        # binding energy
        calc_binding_energy = xps_overrides.pop("calc_binding_energy", False)
        if calc_binding_energy:
            correction_energies = xps_overrides.pop("correction_energies", {})
        # TODO should we override the cutoff_wfc, cutoff_rho by the new pseudo?
        structure_preparation_settings = xps_overrides.pop(
            "structure_preparation_settings", Dict({})
        )
        xps = XpsWorkChain.get_builder_from_protocol(
            code=pw_code,
            structure=structure,
            pseudos=pseudos,
            protocol=protocol,
            elements_list=elements_list,
            calc_binding_energy=Bool(calc_binding_energy),
            correction_energies=Dict(correction_energies),
            core_hole_treatments=core_hole_treatments,
            overrides=xps_overrides,
            structure_preparation_settings=structure_preparation_settings,
            **kwargs,
        )

        xps.pop("relax")
        # xps.pop("structure", None)
        xps.pop("clean_workdir", None)
        builder.xps = xps

        builder.clean_workdir = overrides.get("clean_workdir", Bool(False))
        if "kpoints_distance_override" in overrides:
            builder.kpoints_distance_override = overrides["kpoints_distance_override"]
        if "degauss_override" in overrides:
            builder.degauss_override = overrides["degauss_override"]
        if "smearing_override" in overrides:
            builder.smearing_override = overrides["smearing_override"]

        return builder

    def setup(self):
        """Perform the initial setup of the work chain."""
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_number_of_bands = None
        self.ctx.scf_parent_folder = None

    def should_run_relax(self):
        """Check if the geometry of the input structure should be optimized."""
        return "relax" in self.inputs

    def run_relax(self):
        """Run the `PwRelaxWorkChain`."""
        inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain, namespace="relax"))
        inputs.metadata.call_link_label = "relax"
        inputs.structure = self.ctx.current_structure

        if "kpoints_distance_override" in self.inputs:
            inputs.base.kpoints_distance = self.inputs.kpoints_distance_override
            if "base_scf" in inputs:
                inputs.base_scf.kpoints_distance = self.inputs.kpoints_distance_override

        inputs.base.pw.parameters = inputs.base.pw.parameters.get_dict()
        if "degauss_override" in self.inputs:
            inputs.base.pw.parameters.setdefault("SYSTEM", {})[
                "degauss"
            ] = self.inputs.degauss_override.value

            if "base_scf" in inputs:
                inputs.base_scf_params.pw.parameters = (
                    inputs.base_scf_params.pw.parameters.get_dict()
                )
                inputs.base_scf_params.pw.parameters.setdefault("SYSTEM", {})[
                    "degauss"
                ] = self.inputs.degauss_override.value
        if "smearing_override" in self.inputs:
            inputs.base.pw.parameters.setdefault("SYSTEM", {})[
                "smearing"
            ] = self.inputs.smearing_override.value

            if "base_scf" in inputs:
                inputs.base_scf_params.pw.parameters = (
                    inputs.base_scf_params.pw.parameters.get_dict()
                )
                inputs.base_scf_params.pw.parameters.setdefault("SYSTEM", {})[
                    "smearing"
                ] = self.inputs.smearing_override.value

        inputs = prepare_process_inputs(PwRelaxWorkChain, inputs)
        running = self.submit(PwRelaxWorkChain, **inputs)

        self.report(f"launching PwRelaxWorkChain<{running.pk}>")

        return ToContext(workchain_relax=running)

    def inspect_relax(self):
        """Verify that the `PwRelaxWorkChain` finished successfully."""
        workchain = self.ctx.workchain_relax

        if not workchain.is_finished_ok:
            self.report(
                f"PwRelaxWorkChain failed with exit status {workchain.exit_status}"
            )
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        if "output_structure" in workchain.outputs:
            self.ctx.current_structure = workchain.outputs.output_structure
            self.ctx.current_number_of_bands = (
                workchain.outputs.output_parameters.get_attribute("number_of_bands")
            )
            self.out("structure", self.ctx.current_structure)

    def should_run_bands(self):
        """Check if the band structure should be calculated."""
        return "bands" in self.inputs

    def run_bands(self):
        """Run the `PwBandsWorkChain`."""
        inputs = AttributeDict(self.exposed_inputs(PwBandsWorkChain, namespace="bands"))
        inputs.metadata.call_link_label = "bands"
        inputs.structure = self.ctx.current_structure
        inputs.scf.pw.parameters = inputs.scf.pw.parameters.get_dict()

        if self.ctx.current_number_of_bands:
            inputs.scf.pw.parameters.setdefault("SYSTEM", {}).setdefault(
                "nbnd", self.ctx.current_number_of_bands
            )

        if "kpoints_distance_override" in self.inputs:
            inputs.scf.kpoints_distance = self.inputs.kpoints_distance_override

        if "degauss_override" in self.inputs:
            inputs.scf.pw.parameters.setdefault("SYSTEM", {})[
                "degauss"
            ] = self.inputs.degauss_override.value
        if "smearing_override" in self.inputs:
            inputs.scf.pw.parameters.setdefault("SYSTEM", {})[
                "smearing"
            ] = self.inputs.smearing_override.value

        inputs = prepare_process_inputs(PwBandsWorkChain, inputs)
        running = self.submit(PwBandsWorkChain, **inputs)

        self.report(f"launching PwBandsWorkChain<{running.pk}>")

        return ToContext(workchain_bands=running)

    def inspect_bands(self):
        """Verify that the `PwBandsWorkChain` finished successfully."""
        workchain = self.ctx.workchain_bands

        if not workchain.is_finished_ok:
            self.report(
                f"PwBandsWorkChain failed with exit status {workchain.exit_status}"
            )
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_BANDS

        scf = workchain.get_outgoing(WorkChainNode, link_label_filter="scf").one().node
        self.ctx.scf_parent_folder = scf.outputs.remote_folder
        self.ctx.current_structure = workchain.outputs.primitive_structure

    def should_run_pdos(self):
        """Check if the projected density of states should be calculated."""
        return "pdos" in self.inputs

    def run_pdos(self):
        """Run the `PdosWorkChain`."""
        inputs = AttributeDict(self.exposed_inputs(PdosWorkChain, namespace="pdos"))
        inputs.metadata.call_link_label = "pdos"
        inputs.structure = self.ctx.current_structure
        inputs.nscf.pw.parameters = inputs.nscf.pw.parameters.get_dict()

        if self.ctx.current_number_of_bands:
            inputs.nscf.pw.parameters.setdefault("SYSTEM", {}).setdefault(
                "nbnd", self.ctx.current_number_of_bands
            )

        if self.ctx.scf_parent_folder:
            inputs.pop("scf")
            inputs.nscf.pw.parent_folder = self.ctx.scf_parent_folder
        else:
            if "kpoints_distance_override" in self.inputs:
                inputs.scf.kpoints_distance = self.inputs.kpoints_distance_override

            inputs.scf.pw.parameters = inputs.scf.pw.parameters.get_dict()
            if "degauss_override" in self.inputs:
                inputs.scf.pw.parameters.setdefault("SYSTEM", {})[
                    "degauss"
                ] = self.inputs.degauss_override.value

            if "smearing_override" in self.inputs:
                inputs.scf.pw.parameters.setdefault("SYSTEM", {})[
                    "smearing"
                ] = self.inputs.smearing_override.value

        inputs = prepare_process_inputs(PdosWorkChain, inputs)
        running = self.submit(PdosWorkChain, **inputs)

        self.report(f"launching PdosWorkChain<{running.pk}>")

        return ToContext(workchain_pdos=running)

    def inspect_pdos(self):
        """Verify that the `PdosWorkChain` finished successfully."""
        workchain = self.ctx.workchain_pdos

        if not workchain.is_finished_ok:
            self.report(
                f"PdosWorkChain failed with exit status {workchain.exit_status}"
            )
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PDOS

    def should_run_xspectra(self):
        """Check if the projected density of states should be calculated."""
        self.report(f"should_run_xspectra {'xspectra' in self.inputs}")
        return "xspectra" in self.inputs

    def run_xspectra(self):
        """Run the `XspectraCoreWorkChain`."""
        inputs = AttributeDict(
            self.exposed_inputs(XspectraCoreWorkChain, namespace="xspectra")
        )
        inputs.metadata.call_link_label = "xspectra"
        inputs.structure = self.ctx.current_structure

        # if self.ctx.scf_parent_folder:
        # inputs.pop("scf")
        # inputs.xs_prod.xspectra.parent_folder = self.ctx.scf_parent_folder
        # else:
        #     if "kpoints_distance_override" in self.inputs:
        #         inputs.scf.kpoints_distance = self.inputs.kpoints_distance_override

        #     inputs.scf.pw.parameters = inputs.scf.pw.parameters.get_dict()
        #     if "degauss_override" in self.inputs:
        #         inputs.scf.pw.parameters.setdefault("SYSTEM", {})[
        #             "degauss"
        #         ] = self.inputs.degauss_override.value

        #     if "smearing_override" in self.inputs:
        #         inputs.scf.pw.parameters.setdefault("SYSTEM", {})[
        #             "smearing"
        #         ] = self.inputs.smearing_override.value
        # self.report(f"run_xspectra: {inputs}")
        inputs = prepare_process_inputs(XspectraCoreWorkChain, inputs)
        # self.report(f"run_xspectra: {inputs}")
        running = self.submit(XspectraCoreWorkChain, **inputs)

        self.report(f"launching XspectraCoreWorkChain<{running.pk}>")

        return ToContext(workchain_xspectra=running)

    def inspect_xspectra(self):
        """Verify that the `XspectraCoreWorkChain` finished successfully."""
        workchain = self.ctx.workchain_xspectra

        if not workchain.is_finished_ok:
            self.report(
                f"XspectraCoreWorkChain failed with exit status {workchain.exit_status}"
            )
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_XSPECTRA

    def should_run_xps(self):
        """Check if the band structure should be calculated."""
        return "xps" in self.inputs

    def run_xps(self):
        """Run the `XpsWorkChain`."""
        inputs = AttributeDict(self.exposed_inputs(XpsWorkChain, namespace="xps"))
        inputs.metadata.call_link_label = "xps"
        inputs.structure = self.ctx.current_structure
        inputs.ch_scf.pw.parameters = inputs.ch_scf.pw.parameters.get_dict()

        if self.ctx.current_number_of_bands:
            inputs.ch_scf.pw.parameters.setdefault("SYSTEM", {}).setdefault(
                "nbnd", self.ctx.current_number_of_bands
            )

        if "kpoints_distance_override" in self.inputs:
            inputs.ch_scf.kpoints_distance = self.inputs.kpoints_distance_override

        if "degauss_override" in self.inputs:
            inputs.ch_scf.pw.parameters.setdefault("SYSTEM", {})[
                "degauss"
            ] = self.inputs.degauss_override.value
        if "smearing_override" in self.inputs:
            inputs.ch_scf.pw.parameters.setdefault("SYSTEM", {})[
                "smearing"
            ] = self.inputs.smearing_override.value
        inputs = prepare_process_inputs(XpsWorkChain, inputs)
        running = self.submit(XpsWorkChain, **inputs)

        self.report(f"launching XpsWorkChain<{running.pk}>")

        return ToContext(workchain_xps=running)

    def inspect_xps(self):
        """Verify that the `XpsWorkChain` finished successfully."""
        workchain = self.ctx.workchain_xps

        if not workchain.is_finished_ok:
            self.report(f"XpsWorkChain failed with exit status {workchain.exit_status}")
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_XPS

        # since there are many ch_scf calculation, we don't need to output the folder and structure.
        # scf = workchain.get_outgoing(WorkChainNode, link_label_filter="scf").one().node
        # self.ctx.scf_parent_folder = scf.outputs.remote_folder
        # self.ctx.current_structure = workchain.outputs.primitive_structure

    def results(self):
        """Add the results to the outputs."""
        if "workchain_bands" in self.ctx:
            self.out(
                "band_parameters", self.ctx.workchain_bands.outputs.band_parameters
            )
            self.out("band_structure", self.ctx.workchain_bands.outputs.band_structure)

        if "workchain_pdos" in self.ctx:
            self.out(
                "nscf_parameters",
                self.ctx.workchain_pdos.outputs.nscf.output_parameters,
            )
            self.out("dos", self.ctx.workchain_pdos.outputs.dos.output_dos)
            if "projections_up" in self.ctx.workchain_pdos.outputs.projwfc:
                self.out(
                    "projections_up",
                    self.ctx.workchain_pdos.outputs.projwfc.projections_up,
                )
                self.out(
                    "projections_down",
                    self.ctx.workchain_pdos.outputs.projwfc.projections_down,
                )
            else:
                self.out(
                    "projections", self.ctx.workchain_pdos.outputs.projwfc.projections
                )

        if "workchain_xspectra" in self.ctx:
            self.out(
                "nscf_parameters",
                self.ctx.workchain_xspectra.outputs.nscf.output_parameters,
            )
            self.out("dos", self.ctx.workchain_xspectra.outputs.dos.output_dos)
            if "projections_up" in self.ctx.workchain_xspectra.outputs.projwfc:
                self.out(
                    "projections_up",
                    self.ctx.workchain_xspectra.outputs.projwfc.projections_up,
                )
                self.out(
                    "projections_down",
                    self.ctx.workchain_xspectra.outputs.projwfc.projections_down,
                )
            else:
                self.out(
                    "projections",
                    self.ctx.workchain_xspectra.outputs.projwfc.projections,
                )

        if "workchain_xps" in self.ctx:
            self.out("chemical_shifts", self.ctx.workchain_xps.outputs.chemical_shifts)
            self.out(
                "xps_spectra_cls", self.ctx.workchain_xps.outputs.final_spectra_cls
            )
            self.out(
                "symmetry_analysis_data",
                self.ctx.workchain_xps.outputs.symmetry_analysis_data,
            )
            if "binding_energies" in self.ctx.workchain_xps.outputs:
                self.out(
                    "binding_energies", self.ctx.workchain_xps.outputs.binding_energies
                )
                self.out(
                    "xps_spectra_be", self.ctx.workchain_xps.outputs.final_spectra_be
                )

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report("remote folders will not be cleaned")
            return

        cleaned_calcs = []

        for called_descendant in self.node.called_descendants:
            if isinstance(called_descendant, CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()  # pylint: disable=protected-access
                    cleaned_calcs.append(called_descendant.pk)
                except (OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report(
                f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}"
            )


__version__ = "23.2.0"
