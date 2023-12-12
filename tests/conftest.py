from __future__ import annotations

import io
import pathlib
import tempfile

import pytest
from aiida import orm

pytest_plugins = ["aiida.manage.tests.pytest_fixtures"]


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost


@pytest.fixture
def fixture_code(fixture_localhost):
    """Return an ``InstalledCode`` instance configured to run calculations of given entry point on localhost."""

    def _fixture_code(entry_point_name):
        from aiida.orm import InstalledCode, load_code

        label = f"test.{entry_point_name}"

        try:
            return load_code(label=label)
        except Exception:
            return InstalledCode(
                label=label,
                computer=fixture_localhost,
                filepath_executable="/bin/true",
                default_calc_job_plugin=entry_point_name,
            )

    return _fixture_code


@pytest.fixture
def generate_structure_data():
    """generate a `StructureData` object."""

    def _generate_structure_data(name="silicon", pbc=(True, True, True)):
        if name == "silicon":
            structure = orm.StructureData(
                cell=[
                    [3.84737, 0.0, 0.0],
                    [1.923685, 3.331920, 0.0],
                    [1.923685, 1.110640, 3.141364],
                ]
            )
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Si")
            structure.append_atom(position=(1.923685, 1.110640, 0.785341), symbols="Si")
        elif name == "silica":
            structure = orm.StructureData(
                cell=[
                    [4.18, 0.0, 0.0],
                    [0.0, 4.18, 0.0],
                    [0.0, 0.0, 2.66],
                ]
            )
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Si")
            structure.append_atom(position=(2.09, 2.09, 1.33), symbols="Si")
            structure.append_atom(position=(3.37, 0.81, 1.33), symbols="O")
            structure.append_atom(position=(1.28, 1.28, 0.0), symbols="O")
            structure.append_atom(position=(2.9, 2.9, 0.0), symbols="O")
            structure.append_atom(position=(0.81, 3.37, 1.33), symbols="O")
        structure.pbc = pbc
        return structure

    return _generate_structure_data


@pytest.fixture
def generate_xy_data():
    """Return an ``XyData`` instance."""

    def _generate_xy_data(xvals=None, yvals=None, xlabel=None, ylabel=None):
        """Return an ``XyData`` node.
        xvals and yvals are lists, and should have the same length.
        """
        from aiida.orm import XyData

        xvals = xvals
        yvals = yvals
        xlabel = xlabel
        ylabel = ylabel
        xunits = "n/a"
        yunits = ["n/a"] * len(ylabel)

        xy_node = XyData()
        xy_node.set_x(xvals, xlabel, xunits)
        xy_node.set_y(yvals, ylabel, yunits)
        xy_node.store()
        return xy_node

    return _generate_xy_data


@pytest.fixture
def generate_bands_data():
    """Return a `BandsData` node."""

    def _generate_bands_data():
        """Return a `BandsData` instance with some basic `kpoints` and `bands` arrays."""
        import numpy as np
        from aiida.plugins import DataFactory

        BandsData = DataFactory("core.array.bands")

        kpoints = np.array([[0.0, 0.0, 0.0]])
        bands = np.array([[-5.64024889]])
        bands_data = BandsData()
        bands_data.set_kpoints(kpoints)
        bands_data.set_bands(bands, units="eV")
        bands_data.store()

        return bands_data

    return _generate_bands_data


@pytest.fixture
def generate_projection_data(generate_bands_data):
    """Return an ``ProjectionData`` instance."""

    def _generate_projection_data():
        """Return an ``ProjectionData`` node."""
        import numpy as np
        from aiida.plugins import DataFactory, OrbitalFactory

        ProjectionData = DataFactory("core.array.projection")
        OrbitalCls = OrbitalFactory("core.realhydrogen")

        state_dict = {
            "kind_name": "C",
            "angular_momentum": 0,
            "magnetic_number": 0,
            "radial_nodes": 1,
            "position": [0.0, 0.0, 0.0],
        }
        orbitals = [OrbitalCls(**state_dict)]
        energy_arrays = np.array([1])
        pdos_arrays = np.array([1])

        projection_data = ProjectionData()
        bands_data = generate_bands_data()
        projection_data.set_reference_bandsdata(bands_data)
        projection_data.set_projectiondata(
            orbitals,
            list_of_energy=energy_arrays,
            list_of_pdos=pdos_arrays,
            bands_check=False,
        )
        projection_data.store()
        return projection_data

    return _generate_projection_data


@pytest.fixture(scope="function")
def sssp(aiida_profile, generate_upf_data):
    """Create an SSSP pseudo potential family from scratch."""
    from aiida.common.constants import elements
    from aiida.plugins import GroupFactory

    from aiidalab_qe.common.setup_pseudos import SSSP_VERSION

    aiida_profile.clear_profile()

    SsspFamily = GroupFactory("pseudo.family.sssp")

    cutoffs = {}
    stringency = "standard"

    with tempfile.TemporaryDirectory() as dirpath:
        for values in elements.values():
            element = values["symbol"]

            actinides = (
                "Ac",
                "Th",
                "Pa",
                "U",
                "Np",
                "Pu",
                "Am",
                "Cm",
                "Bk",
                "Cf",
                "Es",
                "Fm",
                "Md",
                "No",
                "Lr",
            )

            if element in actinides:
                continue

            upf = generate_upf_data(element)
            dirpath = pathlib.Path(dirpath)
            filename = dirpath / f"{element}.upf"

            with open(filename, "w+b") as handle:
                with upf.open(mode="rb") as source:
                    handle.write(source.read())
                    handle.flush()

            cutoffs[element] = {
                "cutoff_wfc": 30.0,
                "cutoff_rho": 240.0,
            }

        label = f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"
        family = SsspFamily.create_from_folder(dirpath, label)

    family.set_cutoffs(cutoffs, stringency, unit="Ry")

    return family


@pytest.fixture(scope="function")
def generate_upf_data():
    """Return a `UpfData` instance for the given element a file for which should exist in `tests/fixtures/pseudos`."""

    def _generate_upf_data(element, filename=None):
        """Return `UpfData` node."""
        from aiida_pseudo.data.pseudo import UpfData

        content = f'<UPF version="2.0.1"><PP_HEADER\nelement="{element}"\nz_valence="4.0"\n/></UPF>\n'
        stream = io.BytesIO(content.encode("utf-8"))

        if filename is None:
            filename = f"{element}.upf"

        return UpfData(stream, filename=filename)

    return _generate_upf_data


@pytest.fixture
def pw_code(aiida_local_code_factory):
    """Return a `Code` configured for the pw.x executable."""

    return aiida_local_code_factory(
        label="pw", executable="bash", entry_point="quantumespresso.pw"
    )


@pytest.fixture
def dos_code(aiida_local_code_factory):
    """Return a `Code` configured for the dos.x executable."""
    return aiida_local_code_factory(
        label="dos", executable="bash", entry_point="quantumespresso.dos"
    )


@pytest.fixture
def projwfc_code(aiida_local_code_factory):
    """Return a `Code` configured for the projwfc.x executable."""
    return aiida_local_code_factory(
        label="projwfc",
        executable="bash",
        entry_point="quantumespresso.projwfc",
    )


@pytest.fixture()
def workchain_settings_generator():
    """Return a function that generates a workchain settings dictionary."""
    from aiidalab_qe.app.configuration.workflow import WorkChainSettings

    def _workchain_settings_generator(**kwargs):
        workchain_settings = WorkChainSettings()
        workchain_settings._update_settings(**kwargs)
        return workchain_settings

    return _workchain_settings_generator


@pytest.fixture()
def smearing_settings_generator():
    """Return a function that generates a smearing settings dictionary."""
    from aiidalab_qe.app.configuration.advanced import SmearingSettings

    def _smearing_settings_generator(**kwargs):
        smearing_settings = SmearingSettings()
        smearing_settings.update_settings(**kwargs)
        return smearing_settings

    return _smearing_settings_generator


@pytest.fixture
def app(pw_code, dos_code, projwfc_code):
    from aiidalab_qe.app.main import App

    # Since we use `qe_auto_setup=False`, which will skip the pseudo library installation
    # we need to mock set the installation status to `True` to avoid the blocker message pop up in the
    # submmision step.
    app = App(qe_auto_setup=False)
    app.submit_step.sssp_installation_status.installed = True

    # set up codes
    app.submit_step.pw_code.refresh()
    app.submit_step.codes["dos"].refresh()
    app.submit_step.codes["projwfc"].refresh()

    app.submit_step.pw_code.value = pw_code.uuid
    app.submit_step.codes["dos"].value = dos_code.uuid
    app.submit_step.codes["projwfc"].value = projwfc_code.uuid

    yield app


@pytest.fixture()
@pytest.mark.usefixtures("sssp")
def submit_app_generator(
    app,
    generate_structure_data,
):
    """Return a function that generates a submit step widget."""

    def _submit_app_generator(
        relax_type="positions_cell",
        spin_type="none",
        electronic_type="metal",
        properties=None,
        workchain_protocol="moderate",
        kpoints_distance=0.12,
        smearing="methfessel-paxton",
        degauss=0.015,
        tot_charge=0.0,
        initial_magnetic_moments=0.0,
    ):
        configure_step = app.configure_step
        # Settings
        configure_step.input_structure = generate_structure_data()
        parameters = {
            "workchain": {
                "relax_type": relax_type,
                "spin_type": spin_type,
                "electronic_type": electronic_type,
                "properties": properties or [],
                "protocol": workchain_protocol,
            }
        }
        configure_step.set_configuration_parameters(parameters)
        # Advanced settings
        configure_step.advanced_settings.override.value = True
        configure_step.advanced_settings.total_charge.value = tot_charge
        configure_step.advanced_settings.kpoints_distance.value = kpoints_distance
        configure_step.advanced_settings.magnetization._set_magnetization_values(
            initial_magnetic_moments
        )
        # mimic the behavior of the smearing widget set up
        configure_step.advanced_settings.smearing.smearing.value = smearing
        configure_step.advanced_settings.smearing.degauss.value = degauss
        configure_step.confirm()
        #
        submit_step = app.submit_step
        submit_step.input_structure = generate_structure_data()
        submit_step.resources_config.num_cpus.value = 2

        return app

    return _submit_app_generator


@pytest.fixture
def app_to_submit(app):
    # Step 1: select structure from example
    step1 = app.structure_step
    structure = step1.manager.children[0].children[3]
    structure.children[0].value = structure.children[0].options[1][1]
    step1.confirm()
    # Step 2: configure calculation
    step2 = app.configure_step
    step2.workchain_settings.properties["bands"].run.value = True
    step2.workchain_settings.properties["pdos"].run.value = True
    step2.confirm()
    yield app


@pytest.fixture
def generate_workchain():
    """Generate an instance of a `WorkChain`."""

    def _generate_workchain(process_class, inputs):
        """Generate an instance of a `WorkChain` with the given entry point and inputs.

        :param entry_point: entry point name of the work chain subclass.
        :param inputs: inputs to be passed to process construction.
        :return: a `WorkChain` instance.
        """
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager

        runner = get_manager().get_runner()
        process = instantiate_process(runner, process_class, **inputs)

        return process

    return _generate_workchain


@pytest.fixture
def generate_pdos_workchain(
    fixture_localhost,
    fixture_code,
    generate_xy_data,
    generate_projection_data,
    generate_workchain,
):
    """Generate an instance of a `XpsWorkChain`."""

    def _generate_pdos_workchain(structure, spin_type="none"):
        import numpy as np
        from aiida import engine
        from aiida.orm import Dict, FolderData, RemoteData
        from aiida_quantumespresso.workflows.pdos import PdosWorkChain

        inputs = {
            "pw_code": fixture_code("quantumespresso.pw"),
            "dos_code": fixture_code("quantumespresso.dos"),
            "projwfc_code": fixture_code("quantumespresso.projwfc"),
            "structure": structure,
        }
        builder = PdosWorkChain.get_builder_from_protocol(**inputs)
        inputs = builder._inputs()
        wkchain = generate_workchain(PdosWorkChain, inputs)
        wkchain.setup()
        # wkchain.run_pdos()
        remote = RemoteData(remote_path="/tmp/aiida_run")
        remote.computer = fixture_localhost
        remote.store()
        retrieved = FolderData(tree="/tmp/aiida_run")
        retrieved.store()
        output_parameters = Dict(dict={"fermi_energy": 2.0})
        output_parameters.store()
        proj = generate_projection_data()
        proj.store()
        if spin_type == "none":
            xy = generate_xy_data(
                np.array([1, 2, 3]), [np.array([1, 2, 3])], "X", ["dos"]
            )
            xy.store()
            wkchain.out(
                "dos",
                {
                    "output_dos": xy,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
            wkchain.out(
                "projwfc",
                {
                    "Dos": xy,
                    "projections": proj,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
        else:
            xy = generate_xy_data(
                np.array([1, 2, 3]),
                [np.array([1, 2, 3]), np.array([1, 2, 3])],
                "X",
                ["dos_spin_up", "dos_spin_down"],
            )
            xy.store()
            wkchain.out(
                "dos",
                {
                    "output_dos": xy,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
            wkchain.out(
                "projwfc",
                {
                    "Dos": xy,
                    "projections_up": proj,
                    "projections_down": proj,
                    "output_parameters": output_parameters,
                    "remote_folder": remote,
                    "retrieved": retrieved,
                },
            )
        wkchain.out(
            "nscf",
            {
                "output_parameters": output_parameters,
                "remote_folder": remote,
                "retrieved": retrieved,
            },
        )
        wkchain.update_outputs()
        pdos_node = wkchain.node
        pdos_node.set_exit_status(0)
        pdos_node.set_process_state(engine.ProcessState.FINISHED)
        # set
        return wkchain

    return _generate_pdos_workchain


@pytest.fixture
def generate_bands_workchain(
    fixture_code,
    generate_bands_data,
    generate_workchain,
):
    """Generate an instance of a the WorkChain."""

    def _generate_bands_workchain(structure):
        from copy import deepcopy

        from aiida import engine
        from aiida.orm import Dict
        from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain

        inputs = {
            "code": fixture_code("quantumespresso.pw"),
            "structure": structure,
        }
        builder = PwBandsWorkChain.get_builder_from_protocol(**inputs)
        inputs = builder._inputs()
        inputs["relax"]["base_final_scf"] = deepcopy(inputs["relax"]["base"])
        wkchain = generate_workchain(PwBandsWorkChain, inputs)
        wkchain.setup()
        # run bands and return the process
        output_parameters = Dict(dict={"fermi_energy": 2.0})
        output_parameters.store()
        wkchain.out("scf_parameters", output_parameters)
        wkchain.out("band_parameters", output_parameters)
        #
        band_structure = generate_bands_data()
        band_structure.store()
        wkchain.out("band_structure", band_structure)
        wkchain.update_outputs()
        #
        bands_node = wkchain.node
        bands_node.set_exit_status(0)
        bands_node.set_process_state(engine.ProcessState.FINISHED)
        return wkchain

    return _generate_bands_workchain


@pytest.fixture
def generate_qeapp_workchain(
    app,
    generate_workchain,
    generate_pdos_workchain,
    generate_bands_workchain,
):
    """Generate an instance of the WorkChain."""

    def _generate_qeapp_workchain(
        structure: orm.StructureData | None = None,
        relax_type="positions_cell",
        run_bands=True,
        run_pdos=True,
        spin_type="none",
        initial_magnetic_moments=0.0,
    ):
        from copy import deepcopy

        from aiida.orm.utils.serialize import serialize

        from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
        from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
        from aiidalab_qe.workflows import QeAppWorkChain

        # Step 1: select structure from example
        s1 = app.structure_step
        if structure is None:
            from_example = s1.manager.children[0].children[3]
            # TODO: (unkpcz) using options to set value in test is cranky, instead, use fixture which will make the test more static and robust.
            from_example.children[0].value = from_example.children[0].options[1][1]
        else:
            structure.store()
            aiida_database = s1.manager.children[0].children[2]
            aiida_database.search()
            aiida_database.results.value = structure
        s1.confirm()
        structure = s1.confirmed_structure
        # step 2 configure
        s2: ConfigureQeAppWorkChainStep = app.configure_step
        s2.workchain_settings.relax_type.value = relax_type
        # In order to parepare a complete inputs, I set all the properties to true
        # this can be overrided later
        s2.workchain_settings.properties["bands"].run.value = run_bands
        s2.workchain_settings.properties["pdos"].run.value = run_pdos
        s2.workchain_settings.workchain_protocol.value = "fast"
        s2.workchain_settings.spin_type.value = spin_type
        s2.advanced_settings.magnetization._set_magnetization_values(
            initial_magnetic_moments
        )
        print(s2.advanced_settings.pseudo_family_selector.value)
        s2.confirm()
        # step 3 setup code and resources
        s3: SubmitQeAppWorkChainStep = app.submit_step
        s3.resources_config.num_cpus.value = 4
        builder = s3._create_builder()
        inputs = builder._inputs()
        inputs["relax"]["base_final_scf"] = deepcopy(inputs["relax"]["base"])
        if run_bands:
            inputs["properties"].append("bands")
        if run_pdos:
            inputs["properties"].append("pdos")
        wkchain = generate_workchain(QeAppWorkChain, inputs)
        wkchain.setup()
        # mock output
        if relax_type != "none":
            wkchain.out("structure", s1.confirmed_structure)
        if run_pdos:
            from aiida_quantumespresso.workflows.pdos import PdosWorkChain

            pdos = generate_pdos_workchain(structure, spin_type)
            wkchain.out_many(
                wkchain.exposed_outputs(pdos.node, PdosWorkChain, namespace="pdos")
            )
        if run_bands:
            from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain

            bands = generate_bands_workchain(structure)
            wkchain.out_many(
                wkchain.exposed_outputs(bands.node, PwBandsWorkChain, namespace="bands")
            )
        wkchain.update_outputs()
        # set ui_parameters
        qeapp_node = wkchain.node
        qeapp_node.base.extras.set("ui_parameters", serialize(s3.ui_parameters))
        return wkchain

    return _generate_qeapp_workchain
