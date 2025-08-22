import pytest

from aiida import orm
from aiidalab_qe.app.configuration.advanced.pseudos import (
    PseudosConfigurationSettingsModel,
    PseudosConfigurationSettingsPanel,
)
from aiidalab_qe.app.configuration.advanced.pseudos.uploader import (
    PseudoPotentialUploader,
    PseudoPotentialUploaderModel,
)
from aiidalab_qe.setup.pseudos import (
    PSEUDODOJO_VERSION,
    SSSP_VERSION,
    _construct_cmd,
    _install_pseudos,
    pseudos_to_install,
)


def test_setup_pseudos_cmd(tmp_path):
    """Test _construct_cmd function in setup.pseudos"""

    # SSSP family
    pseudo_family = f"SSSP/{SSSP_VERSION}/PBE/efficiency"
    cmd = _construct_cmd(pseudo_family)
    assert cmd == [
        "aiida-pseudo",
        "install",
        "sssp",
        "--functional",
        "PBE",
        "--version",
        f"{SSSP_VERSION}",
        "-p",
        "efficiency",
    ]

    # PseudoDojo family
    pseudo_family = f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf"
    cmd = _construct_cmd(pseudo_family)
    assert cmd == [
        "aiida-pseudo",
        "install",
        "pseudo-dojo",
        "--functional",
        "PBEsol",
        "--version",
        f"{PSEUDODOJO_VERSION}",
        "-p",
        "standard",
        "--relativistic",
        "SR",
        "--pseudo-format",
        "upf",
    ]

    # with download_only option
    pseudo_family = f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf"
    cmd = _construct_cmd(pseudo_family, download_only=True)
    assert cmd == [
        "aiida-pseudo",
        "install",
        "pseudo-dojo",
        "--functional",
        "PBEsol",
        "--version",
        f"{PSEUDODOJO_VERSION}",
        "-p",
        "standard",
        "--relativistic",
        "SR",
        "--pseudo-format",
        "upf",
        "--download-only",
    ]

    # with cwd option to specify the source folder
    pseudo_family = f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf"
    cmd = _construct_cmd(pseudo_family, cwd=tmp_path)

    # since the source file not exist, the cmd should be the same as above
    assert "--from-download" not in cmd

    # mock the source file
    source_file = (
        tmp_path
        / f"PseudoDojo_{PSEUDODOJO_VERSION}_PBEsol_SR_standard_upf.aiida_pseudo"
    )
    source_file.touch()
    cmd = _construct_cmd(pseudo_family, cwd=tmp_path)
    assert cmd == [
        "aiida-pseudo",
        "install",
        "pseudo-dojo",
        "--functional",
        "PBEsol",
        "--version",
        f"{PSEUDODOJO_VERSION}",
        "-p",
        "standard",
        "--relativistic",
        "SR",
        "--pseudo-format",
        "upf",
        "--from-download",
        f"{tmp_path!s}/PseudoDojo_{PSEUDODOJO_VERSION}_PBEsol_SR_standard_upf.aiida_pseudo",
    ]


@pytest.mark.slow
@pytest.mark.usefixtures("aiida_profile_clean")
def test_pseudos_installation():
    """Test install_pseudos"""
    # Test by compare the pseudos_to_install before and after the installation
    assert len(pseudos_to_install()) == 12
    EXPECTED_PSEUDOS = {
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/standard/upf",
        f"SSSP/{SSSP_VERSION}/PBE/efficiency",
    }

    # Install the pseudos
    list(_install_pseudos(EXPECTED_PSEUDOS))

    # Two pseudos are installed
    assert len(pseudos_to_install()) == 10


@pytest.mark.slow
@pytest.mark.usefixtures("aiida_profile_clean")
def test_download_and_install_pseudo_from_file(tmp_path):
    """Test download and install pseudo from file."""
    assert len(pseudos_to_install()) == 12
    EXPECTED_PSEUDOS = {
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/standard/upf",
        f"SSSP/{SSSP_VERSION}/PBE/efficiency",
    }

    # Download the pseudos to the tmp_path but not install
    list(_install_pseudos(EXPECTED_PSEUDOS, download_only=True, cwd=tmp_path))

    assert len(pseudos_to_install()) == 12
    assert len(list(tmp_path.iterdir())) == 2

    # Install the pseudos from the tmp_path
    list(_install_pseudos(EXPECTED_PSEUDOS, cwd=tmp_path))

    # Two pseudos are installed
    assert len(pseudos_to_install()) == 10


def test_pseudos_settings(generate_structure_data, generate_upf_data):
    model = PseudosConfigurationSettingsModel()
    pseudos = PseudosConfigurationSettingsPanel(model=model)

    silicon = generate_structure_data("silicon")
    model.input_structure = silicon

    # Test the default family
    assert model.family == f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"
    assert "Si" in model.dictionary.keys()
    assert model.ecutwfc == 30
    assert model.ecutrho == 240

    # Test protocol-dependent family change
    model.protocol = "stringent"
    assert model.family == f"SSSP/{SSSP_VERSION}/PBEsol/precision"

    # Test functional-dependent family change
    model.functional = "PBE"
    assert model.family == f"SSSP/{SSSP_VERSION}/PBE/precision"

    # Test library-dependent family change
    model.library = "PseudoDojo stringent"
    assert model.family == f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/stringent/upf"

    # Test spin-orbit-dependent family change
    model.spin_orbit = "soc"
    assert model.family == f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/FR/stringent/upf"

    # Reset the external dependencies of the model
    model.spin_orbit = "wo_soc"
    model.protocol = "balanced"

    # Test that changing the structure triggers a reset
    silica = generate_structure_data("silica")
    model.input_structure = silica
    assert model.family == f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"
    assert "Si" in model.dictionary.keys()
    assert "O" in model.dictionary.keys()

    pseudos.render()

    # Check uploaders
    assert len(pseudos.setter_widget.children) == 2

    message = "{ecutwfc} | {ecutrho} | {functional} | {relativistic}"

    # Check Si uploader (Si.upf)
    Si_uploader: PseudoPotentialUploader = pseudos.setter_widget.children[0]
    assert Si_uploader._model.kind_name == "Si"
    assert Si_uploader._model.kind_symbol == "Si"
    pseudo = orm.load_node(model.dictionary["Si"])
    assert Si_uploader._model.pseudo == pseudo
    assert Si_uploader._model.cutoffs == [30, 240]
    assert Si_uploader.pseudo_filename.value == pseudo.filename
    assert (
        message.format(
            ecutwfc=30.0,
            ecutrho=240.0,
            functional="PBEsol",
            relativistic="N/A",
        )
        in Si_uploader._model.info
    )

    # Check O uploader (O.upf)
    O_uploader: PseudoPotentialUploader = pseudos.setter_widget.children[1]
    assert O_uploader._model.kind_name == "O"
    assert O_uploader._model.kind_symbol == "O"
    pseudo = orm.load_node(model.dictionary["O"])
    assert O_uploader._model.pseudo == pseudo
    assert O_uploader._model.cutoffs == [30, 240]
    assert O_uploader.pseudo_filename.value == pseudo.filename
    assert (
        message.format(
            ecutwfc=30.0,
            ecutrho=240.0,
            functional="PBEsol",
            relativistic="N/A",
        )
        in O_uploader._model.info
    )

    # Test reset from uploaded state
    new_O_pseudo = generate_upf_data(
        "O",
        "O_new.upf",
        params={"accuracy": "high"},
    )
    O_uploader._on_file_upload(
        {
            "new": {
                "O_new.upf": {
                    "content": bytes(
                        new_O_pseudo.get_content(),
                        encoding="utf-8",
                    ),
                },
            },
        }
    )
    assert O_uploader._model.pseudo.filename == "O_new.upf"
    assert O_uploader.pseudo_filename.value == "O_new.upf"

    model.reset()

    assert O_uploader._model.pseudo.filename != "O_new.upf"
    assert O_uploader.pseudo_filename.value != "O_new.upf"


def test_pseudo_upload_widget(generate_upf_data):
    """Test the pseudo upload widget."""

    # Test that the kind can be not the element symbol
    # the widget initialize with the pseudo as input to mock how it will
    # be used in PseudoSetter when the pseudo family is set.
    old_pseudo = generate_upf_data("O", "O.upf")
    old_pseudo.store()

    model = PseudoPotentialUploaderModel(kind_name="O1", kind_symbol="O")
    uploader = PseudoPotentialUploader(model=model)

    model.pseudo = old_pseudo
    model.cutoffs = [30, 240]
    model.update_pseudo_info()

    message = "{ecutwfc} | {ecutrho}"

    assert model.pseudo.filename == "O.upf"
    assert message.format(ecutwfc=30.0, ecutrho=240.0) in model.info
    assert not model.message

    uploader.render()

    # Check different element is rejected
    different_element = generate_upf_data("Si", "Si.upf")
    uploader._on_file_upload(
        {
            "new": {
                "Si.upf": {
                    "content": bytes(
                        different_element.get_content(),
                        encoding="utf-8",
                    )
                },
            },
        }
    )
    assert model.pseudo.filename == "O.upf"
    assert "does not match" in model.message

    # Check identical content is rejected in favor of existing one
    same_content_different_name = generate_upf_data("O", "O_copy.upf")
    uploader._on_file_upload(
        {
            "new": {
                "O_copy.upf": {
                    "content": bytes(
                        same_content_different_name.get_content(),
                        encoding="utf-8",
                    ),
                },
            },
        }
    )
    assert model.pseudo.filename == "O.upf"
    assert "Identical pseudo" in model.message

    # Check different content but same filename is rejected
    different_content_same_filename = generate_upf_data(
        "O",
        "O.upf",
        params={"accuracy": "low"},
    )
    uploader._on_file_upload(
        {
            "new": {
                "O.upf": {
                    "content": bytes(
                        different_content_same_filename.get_content(),
                        encoding="utf-8",
                    ),
                },
            },
        }
    )
    assert model.pseudo.filename == "O.upf"
    assert "rename your file" in model.message

    # Check invalid pseudo content is rejected
    uploader._on_file_upload(
        {
            "new": {
                "O_invalid.upf": {
                    "content": b"<UPF version='...'>...</UPF>",
                },
            },
        }
    )
    assert "not a valid UPF file" in model.message

    # Check valid pseudo is accepted
    valid = generate_upf_data(
        "O",
        "O_valid.upf",
        params={"accuracy": "normal"},
    )
    uploader._on_file_upload(
        {
            "new": {
                "O_valid.upf": {
                    "content": bytes(
                        valid.get_content(),
                        encoding="utf-8",
                    ),
                },
            },
        }
    )
    assert model.pseudo.filename == "O_valid.upf"
    assert "uploaded successfully" in model.message


def test_missing_pseudos(generate_structure_data):
    """Test that the model handles missing pseudos correctly."""
    model = PseudosConfigurationSettingsModel()
    _ = PseudosConfigurationSettingsPanel(model=model)
    model.input_structure = generate_structure_data("CeO")
    model.functional = "PBEsol"
    model.library = "PseudoDojo standard"
    assert model.family == "PseudoDojo/0.4/PBEsol/SR/standard/upf"
    assert model.dictionary["Ce"] is None
    assert model.dictionary["O"] is not None
    assert len(model.blockers) == 1
    assert "does not contain a pseudopotential for Ce" in model.blockers[0]


def test_functional_mismatch_blocker(generate_structure_data):
    """Test blocker for inconsistent functional across selected pseudopotentials."""
    model = PseudosConfigurationSettingsModel()
    _ = PseudosConfigurationSettingsPanel(model=model)
    model.input_structure = generate_structure_data("silica")
    model.functionals = ["PBE", "PBEsol"]
    assert len(model.blockers) == 1
    assert "must have the same exchange-correlation" in model.blockers[0]


def test_relativistic_mismatch_blocker(generate_structure_data):
    """Test blocker for inconsistent relativistic treatment across selected
    pseudopotentials with and without SOC.
    """
    model = PseudosConfigurationSettingsModel()
    _ = PseudosConfigurationSettingsPanel(model=model)
    model.input_structure = generate_structure_data("silica")

    # Check with SOC
    model.spin_orbit = "soc"
    model.family = "SSSP/1.3/PBEsol/efficiency"
    assert len(model.blockers) == 1
    assert "pseudopotentials must be fully relativistic" in model.blockers[0]
