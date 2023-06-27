# AiiDA imports.
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, if_
from aiida.plugins import DataFactory

# AiiDA Quantum ESPRESSO plugin inputs.
from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain

XyData = DataFactory("core.array.xy")
StructureData = DataFactory("core.structure")
BandsData = DataFactory("core.array.bands")
Orbital = DataFactory("core.orbital")


class QeAppWorkChain(WorkChain):
    """WorkChain designed to calculate the requested properties in the AiiDAlab Quantum ESPRESSO app."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('structure', valid_type=StructureData,
                   help='The inputs structure.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False),
                   help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.expose_inputs(PwRelaxWorkChain, namespace='relax', exclude=('clean_workdir', 'structure', 'base_final_scf'),
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
            cls.results
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
                       message='The PwRelaxWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_BANDS',
                       message='The PwBandsWorkChain sub process failed')
        spec.exit_code(404, 'ERROR_SUB_PROCESS_FAILED_PDOS',
                       message='The PdosWorkChain sub process failed')
        spec.output('structure', valid_type=StructureData, required=False)
        spec.output('band_parameters', valid_type=orm.Dict, required=False)
        spec.output('band_structure', valid_type=BandsData, required=False)
        spec.output('nscf_parameters', valid_type=orm.Dict, required=False)
        spec.output('dos', valid_type=XyData, required=False)
        spec.output('projections', valid_type=Orbital, required=False)
        spec.output('projections_up', valid_type=Orbital, required=False)
        spec.output('projections_down', valid_type=Orbital, required=False)
        # yapf: enable

    @classmethod
    def get_builder_from_protocol(cls, structure, parameters):
        """Return a builder prepopulated with inputs selected according to the chosen protocol."""
        from aiidalab_qe.app.plugins.bands.workchain import (
            get_builder as get_bands_builder,
        )
        from aiidalab_qe.app.plugins.pdos.workchain import (
            get_builder as get_pdos_builder,
        )
        from aiidalab_qe.app.plugins.relax.workchain import (
            get_builder as get_relax_builder,
        )

        builder = cls.get_builder()
        # Set the structure.
        builder.structure = structure
        codes = parameters.pop("codes", {})
        #
        # TODO do we always need to run a relax workchain
        relax_builder = get_relax_builder(codes, structure, parameters)
        builder.relax = relax_builder

        # bands workchain
        if parameters["workflow"]["properties"]["bands"]:
            bands_builder = get_bands_builder(codes, structure, parameters)
            builder.bands = bands_builder

        # pdos workchain
        if parameters["workflow"]["properties"]["pdos"]:
            pdos_builder = get_pdos_builder(codes, structure, parameters)
            builder.pdos = pdos_builder

        # TODO check if we need to clean the workdir
        # builder.clean_workdir = orm.Bool(clean_workdir)

        return builder

    def setup(self):
        """Perform the initial setup of the work chain, setup the input structure
        and the logic to determine which sub workchains to run.
        """
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

        scf = (
            workchain.get_outgoing(orm.WorkChainNode, link_label_filter="scf")
            .one()
            .node
        )
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

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report("remote folders will not be cleaned")
            return

        cleaned_calcs = []

        for called_descendant in self.node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()  # pylint: disable=protected-access
                    cleaned_calcs.append(called_descendant.pk)
                except (OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report(
                f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}"
            )
