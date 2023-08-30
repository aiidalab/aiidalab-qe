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
def generate_structure_data():
    """generate a `StructureData` object."""

    def _generate_structure_data(name="silicon"):
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

        return structure

    return _generate_structure_data


@pytest.fixture(scope="session", autouse=True)
def sssp(aiida_profile, generate_upf_data):
    """Create an SSSP pseudo potential family from scratch."""
    from aiida.common.constants import elements
    from aiida.plugins import GroupFactory

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

        label = "SSSP/1.2/PBEsol/efficiency"
        family = SsspFamily.create_from_folder(dirpath, label)

    family.set_cutoffs(cutoffs, stringency, unit="Ry")

    return family


@pytest.fixture(scope="session")
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
        label="pw@localhost", executable="bash", entry_point="quantumespresso.pw"
    )


@pytest.fixture
def dos_code(aiida_local_code_factory):
    """Return a `Code` configured for the dos.x executable."""
    return aiida_local_code_factory(
        label="dos@localhost", executable="bash", entry_point="quantumespresso.dos"
    )


@pytest.fixture
def projwfc_code(aiida_local_code_factory):
    """Return a `Code` configured for the projwfc.x executable."""
    return aiida_local_code_factory(
        label="projwfc@localhost",
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
def initial_magnetic_moments_generator(generate_structure_data):
    """Retturn a function that generatates a initial_magnetic_moments dictionary"""
    from aiidalab_qe.app.configuration.advanced import MagnetizationSettings

    def _initial_moments_generator(**kwargs):
        initial_magnetic_moments = MagnetizationSettings()
        initial_magnetic_moments.input_structure = generate_structure_data()
        initial_magnetic_moments.update_kinds_widget()
        initial_magnetic_moments._set_magnetization_values(**kwargs)
        return initial_magnetic_moments

    return _initial_moments_generator


@pytest.fixture()
def tot_charge_generator():
    """Return a function that generates a tot_charge dictionary."""
    from aiidalab_qe.app.configuration.advanced import TotalCharge

    def _tot_charge_generator(**kwargs):
        tot_charge = TotalCharge()
        tot_charge._update_settings(**kwargs)
        return tot_charge

    return _tot_charge_generator


@pytest.fixture()
def smearing_settings_generator():
    """Return a function that generates a smearing settings dictionary."""
    from aiidalab_qe.app.configuration.advanced import SmearingSettings

    def _smearing_settings_generator(**kwargs):
        smearing_settings = SmearingSettings()
        smearing_settings._update_settings(**kwargs)
        return smearing_settings

    return _smearing_settings_generator


@pytest.fixture()
def kpoints_settings_generator():
    """Return a function that generates a kpoints settings dictionary."""
    from aiidalab_qe.app.configuration.advanced import KpointSettings

    def _kpoints_settings_generator(**kwargs):
        kpoints_settings = KpointSettings()
        kpoints_settings._update_settings(**kwargs)
        return kpoints_settings

    return _kpoints_settings_generator


@pytest.fixture()
@pytest.mark.usefixtures("sssp")
def submit_step_widget_generator(
    pw_code,
    dos_code,
    projwfc_code,
    generate_structure_data,
    workchain_settings_generator,
    smearing_settings_generator,
    kpoints_settings_generator,
    tot_charge_generator,
    initial_magnetic_moments_generator,
):
    """Return a function that generates a submit step widget."""
    from aiidalab_qe.app.configuration.advanced import AdvancedSettings
    from aiidalab_qe.app.configuration.pseudos import PseudoFamilySelector, PseudoSetter
    from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep

    def _submit_step_widget_generator(
        relax_type="positions_cell",
        spin_type="none",
        electronic_type="metal",
        bands_run=True,
        pdo_run=True,
        workchain_protocol="moderate",
        kpoints_distance=0.12,
        smearing="methfessel-paxton",
        degauss=0.015,
        override_protocol_smearing=True,
        tot_charge=0.0,
        initial_magnetic_moments=0.0,
    ):
        submit_step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
        submit_step.input_structure = generate_structure_data()
        submit_step.pseudo_family_selector = PseudoFamilySelector()
        submit_step.pseudo_setter = PseudoSetter()

        submit_step.pw_code.value = pw_code.uuid
        submit_step.dos_code.value = dos_code.uuid
        submit_step.projwfc_code.value = projwfc_code.uuid

        # Settings
        submit_step.workchain_settings = workchain_settings_generator(
            relax_type=relax_type,
            spin_type=spin_type,
            electronic_type=electronic_type,
            bands_run=bands_run,
            pdos_run=pdo_run,
            workchain_protocol=workchain_protocol,
        )
        # Advanced settings
        submit_step.advanced_settings = AdvancedSettings()
        submit_step.advanced_settings.override.value = True

        submit_step.advanced_settings.tot_charge = tot_charge_generator(
            tot_charge=tot_charge,
        )

        submit_step.advanced_settings.magnetization = (
            initial_magnetic_moments_generator(
                initial_magnetic_moments=initial_magnetic_moments
            )
        )

        submit_step.advanced_settings.kpoints = kpoints_settings_generator(
            kpoints_distance=kpoints_distance
        )
        submit_step.advanced_settings.smearing = smearing_settings_generator(
            smearing=smearing,
            degauss=degauss,
        )
        submit_step.advanced_settings.smearing.override.value = (
            override_protocol_smearing
        )

        return submit_step

    return _submit_step_widget_generator
