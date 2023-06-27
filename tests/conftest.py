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


@pytest.fixture()
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


@pytest.fixture()
@pytest.mark.usefixtures("sssp")
def submit_step_widget_generator(
    pw_code,
    dos_code,
    projwfc_code,
    structure_data_object,
):
    """Return a function that generates a submit step widget."""
    from aiidalab_qe.app.submit import SubmitQeAppWorkChainStep

    # I removed all the parameters related with configure step.
    def _submit_step_widget_generator():
        submit_step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
        submit_step.input_structure = structure_data_object

        submit_step.pw_code.value = pw_code.uuid
        submit_step.dos_code.value = dos_code.uuid
        submit_step.projwfc_code.value = projwfc_code.uuid

        return submit_step

    return _submit_step_widget_generator


# I try to use the usefixtures decorator but it does not work
# so I pass the pw_code etc to the parameters list
@pytest.fixture
def app(pw_code, dos_code, projwfc_code, sssp):
    from aiidalab_qe.app.app import QEApp

    app = QEApp(qe_auto_setup=False)
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
def qeapp_workchain(app, generate_workchain):
    """Generate an instance of a `XpsWorkChain`."""
    from aiida import engine

    from aiidalab_qe.workflows import QeAppWorkChain

    # Step 1: select structure from example
    s1 = app.steps.steps[0][1]
    structure = s1.manager.children[0].children[3]
    structure.children[0].value = structure.children[0].options[1][1]
    s1.confirm()
    # step 2
    s2 = app.steps.steps[1][1]
    s2.workchain_settings.relax_type.value = "positions_cell"
    s2.workchain_settings.properties["bands"].run.value = True
    s2.workchain_settings.properties["pdos"].run.value = True
    s2.basic_settings.workchain_protocol.value = "fast"
    s2.confirm()
    # step 3
    #
    s3 = app.steps.steps[2][1]
    s3.resources_config.num_cpus.value = 4
    builder, ui_parameters = s3._create_builder()

    qeapp_process = generate_workchain(QeAppWorkChain, builder._inputs())
    qeapp_node = qeapp_process.node
    qeapp_node.set_exit_status(0)
    qeapp_node.set_process_state(engine.ProcessState.FINISHED)
    qeapp_node.store()
    # set
    qeapp_node.base.extras.set("ui_parameters", ui_parameters)

    return qeapp_process
