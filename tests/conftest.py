import io
import pathlib
import tempfile

import pytest
from aiida import plugins

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
def structure_data_object():
    """Return a `StructureData` object."""
    StructureData = plugins.DataFactory("core.structure")  # noqa: N806
    structure = StructureData(
        cell=[
            [3.84737, 0.0, 0.0],
            [1.923685, 3.331920, 0.0],
            [1.923685, 1.110640, 3.141364],
        ]
    )
    structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Si")
    structure.append_atom(position=(1.923685, 1.110640, 0.785341), symbols="Si")

    return structure


@pytest.fixture
def generate_xy_data():
    """Return an ``XyData`` instance."""

    def _generate_xy_data():
        """Return an ``XyData`` node."""
        import numpy as np
        from aiida.orm import XyData

        xvals = [1, 2, 3]
        yvals = [10, 20, 30]
        xlabel = "X"
        ylabel = "dos"
        xunits = "n/a"
        yunits = "n/a"

        xy_node = XyData()
        xy_node.set_x(np.array(xvals), xlabel, xunits)
        xy_node.set_y([np.array(yvals)], [ylabel], [yunits])
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

        kpoints = np.array([[0.0, 0.0, 0.0]])
        bands = np.array([[-5.64024889]])
        BandsData = DataFactory("core.array.bands")
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
        from aiida.orm import ProjectionData
        from aiida.plugins import OrbitalFactory

        OrbitalCls = OrbitalFactory("core.realhydrogen")
        state_dict = {
            "kind_name": "C",
            "angular_momentum": 0,
            "magnetic_number": 0,
            "radial_nodes": 1,
            "position": [0.0, 0.0, 0.0],
        }
        orbitals = [OrbitalCls(**state_dict)]
        # projections = np.array([[1]])
        energy_arrays = np.array([1])
        pdos_arrays = np.array([1])

        projection_data = ProjectionData()
        bands_data = generate_bands_data()
        projection_data.set_reference_bandsdata(bands_data)
        projection_data.set_projectiondata(
            orbitals,
            # list_of_projections=projections,
            list_of_energy=energy_arrays,
            list_of_pdos=pdos_arrays,
            bands_check=False,
        )
        projection_data.store()
        return projection_data

    return _generate_projection_data


@pytest.fixture(scope="session", autouse=True)
def sssp(aiida_profile, generate_upf_data):
    """Create an SSSP pseudo potential family from scratch."""
    from aiida.plugins import GroupFactory
    from aiida_pseudo.data.pseudo import UpfData

    from aiidalab_qe.app.sssp import EXPECTED_PSEUDOS

    aiida_profile.clear_profile()

    SsspFamily = GroupFactory("pseudo.family.sssp")
    PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")

    cutoffs = {}
    stringency = "standard"

    with tempfile.TemporaryDirectory() as dirpath:
        for element in ["Si"]:
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
        for label in EXPECTED_PSEUDOS:
            if label.startswith("SSSP"):
                family = SsspFamily.create_from_folder(dirpath, label)
            else:
                family = PseudoDojoFamily.create_from_folder(
                    dirpath, label, pseudo_type=UpfData
                )
            family.set_cutoffs(cutoffs, stringency, unit="Ry")

    return family


@pytest.fixture(scope="session")
def generate_upf_data():
    """Return a `UpfData` instance for the given element a file for which should exist in `tests/fixtures/pseudos`."""

    def _generate_upf_data(element):
        """Return `UpfData` node."""
        from aiida_pseudo.data.pseudo import UpfData

        content = f'<UPF version="2.0.1"><PP_HEADER\nelement="{element}"\nz_valence="4.0"\n/></UPF>\n'
        stream = io.BytesIO(content.encode("utf-8"))
        return UpfData(stream, filename=f"{element}.upf")

    return _generate_upf_data


@pytest.fixture
def pw_code(aiida_local_code_factory):
    """Return a `Code` configured for the pw.x executable."""
    from aiida.common import exceptions
    from aiida.orm import load_code

    try:
        return load_code(label="pw-7.2")
    except exceptions.NotExistent:
        return aiida_local_code_factory(
            label="pw-7.2", executable="bash", entry_point="quantumespresso.pw"
        )


@pytest.fixture
def pw_code_70(aiida_local_code_factory):
    """Return a `Code` configured for the pw.x executable."""
    return aiida_local_code_factory(
        label="pw-7.0", executable="bash", entry_point="quantumespresso.pw"
    )


@pytest.fixture
def dos_code(aiida_local_code_factory):
    """Return a `Code` configured for the dos.x executable."""
    return aiida_local_code_factory(
        label="dos-7.2", executable="bash", entry_point="quantumespresso.dos"
    )


@pytest.fixture
def projwfc_code(aiida_local_code_factory):
    """Return a `Code` configured for the projwfc.x executable."""
    return aiida_local_code_factory(
        label="projwfc-7.2",
        executable="bash",
        entry_point="quantumespresso.projwfc",
    )


@pytest.fixture
def workchain_settings_generator():
    """Return a function that generates a workchain settings dictionary."""
    from aiidalab_qe.app.configure.workflow import WorkChainSettings

    def _workchain_settings_generator(**kwargs):
        workchain_settings = WorkChainSettings()
        workchain_settings._update_settings(**kwargs)
        return workchain_settings

    return _workchain_settings_generator


# I removed all the fixture for kinpoints, smearing, tot_charge
# because there are not a seperate class anymore. They are just
# a part of AdvancedSettings class.


# I try to use the usefixtures decorator but it does not work
# so I pass the pw_code etc to the parameters list
@pytest.fixture
def app(pw_code, dos_code, projwfc_code, sssp):
    from aiidalab_qe.app.app import QEApp

    app = QEApp(qe_auto_setup=False)

    yield app


@pytest.fixture
def app_to_submit(app):
    # Step 1: select structure from example
    step1 = app.steps.steps[0][1]
    structure = step1.manager.children[0].children[3]
    structure.children[0].value = structure.children[0].options[1][1]
    step1.confirm()
    # Step 2: configure calculation
    step2 = app.steps.steps[1][1]
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
def generate_qeapp_workchain(app, generate_workchain):
    """Generate an instance of a `XpsWorkChain`."""

    def _generate_qeapp_workchain(
        relax_type="positions_cell", run_bands=True, run_pdos=True
    ):
        from aiida import engine

        from aiidalab_qe.workflows import QeAppWorkChain

        # Step 1: select structure from example
        s1 = app.steps.steps[0][1]
        structure = s1.manager.children[0].children[3]
        structure.children[0].value = structure.children[0].options[1][1]
        s1.confirm()
        # step 2 configure
        s2 = app.steps.steps[1][1]
        s2.workchain_settings.relax_type.value = relax_type
        # In order to parepare a complete inputs, I set all the properties to true
        # I wil override this later
        s2.workchain_settings.properties["bands"].run.value = True
        s2.workchain_settings.properties["pdos"].run.value = True
        s2.basic_settings.workchain_protocol.value = "fast"
        s2.confirm()
        # step 3 setup code and resources
        #
        s3 = app.steps.steps[2][1]
        s3.resources_config.num_cpus.value = 4
        builder, ui_parameters = s3._create_builder()
        inputs = builder._inputs()
        # override the workflow
        inputs["properties"]["bands"] = run_bands
        inputs["properties"]["pdos"] = run_pdos
        qeapp_process = generate_workchain(QeAppWorkChain, inputs)
        qeapp_node = qeapp_process.node
        qeapp_node.set_exit_status(0)
        qeapp_node.set_process_state(engine.ProcessState.FINISHED)
        # set
        qeapp_node.base.extras.set("ui_parameters", ui_parameters)
        return qeapp_process

    return _generate_qeapp_workchain


@pytest.fixture
def generate_pdos_workchain(
    structure_data_object,
    fixture_localhost,
    fixture_code,
    generate_xy_data,
    generate_projection_data,
    generate_workchain,
):
    """Generate an instance of a `XpsWorkChain`."""

    def _generate_pdos_workchain():
        from aiida import engine
        from aiida.orm import Dict, FolderData, RemoteData
        from aiida_quantumespresso.workflows.pdos import PdosWorkChain

        inputs = {
            "pw_code": fixture_code("quantumespresso.pw"),
            "dos_code": fixture_code("quantumespresso.dos"),
            "projwfc_code": fixture_code("quantumespresso.projwfc"),
            "structure": structure_data_object,
        }
        builder = PdosWorkChain.get_builder_from_protocol(**inputs)
        inputs = builder._inputs()
        wkchain = generate_workchain(PdosWorkChain, inputs)
        wkchain.setup()
        # run pdos and return the process
        xy = generate_xy_data()
        xy.store()
        remote = RemoteData(remote_path="/tmp/aiida_run")
        remote.computer = fixture_localhost
        remote.store()
        retrieved = FolderData(tree="/tmp/aiida_run")
        retrieved.store()
        output_parameters = Dict(dict={"fermi_energy": 2.0})
        output_parameters.store()
        wkchain.out(
            "dos",
            {
                "output_dos": xy,
                "output_parameters": output_parameters,
                "remote_folder": remote,
                "retrieved": retrieved,
            },
        )
        proj = generate_projection_data()
        proj.store()
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
    structure_data_object,
    fixture_localhost,
    fixture_code,
    generate_xy_data,
    generate_bands_data,
    generate_workchain,
):
    """Generate an instance of a `XpsWorkChain`."""

    def _generate_bands_workchain():
        from copy import deepcopy

        from aiida import engine
        from aiida.orm import Dict
        from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain

        inputs = {
            "code": fixture_code("quantumespresso.pw"),
            "structure": structure_data_object,
        }
        builder = PwBandsWorkChain.get_builder_from_protocol(**inputs)
        inputs = builder._inputs()
        inputs["relax"]["base_final_scf"] = deepcopy(inputs["relax"]["base"])
        wkchain = generate_workchain(PwBandsWorkChain, inputs)
        wkchain.setup()
        # run pdos and return the process
        output_parameters = Dict(dict={"fermi_energy": 2.0})
        output_parameters.store()
        wkchain.out("scf_parameters", output_parameters)
        wkchain.out("band_parameters", output_parameters)
        #
        band_structure = generate_bands_data()
        band_structure.store()
        wkchain.out("band_structure", band_structure)
        wkchain.update_outputs()
        bands_node = wkchain.node
        bands_node.set_exit_status(0)
        bands_node.set_process_state(engine.ProcessState.FINISHED)
        # set
        return wkchain

    return _generate_bands_workchain
