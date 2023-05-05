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
