import pytest

from aiida import orm
from aiidalab_qe.app.configuration.advanced.pseudos import PseudoUploadWidget
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
    from aiidalab_qe.app.configuration.advanced.pseudos import (
        PseudosConfigurationSettingsModel,
        PseudosConfigurationSettingsPanel,
    )

    model = PseudosConfigurationSettingsModel()
    pseudos = PseudosConfigurationSettingsPanel(model=model)

    # Test the default family
    model.spin_orbit = "wo_soc"
    assert model.family == f"SSSP/{SSSP_VERSION}/PBEsol/efficiency"

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
    model.protocol = "balanced"
    assert model.family == f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/FR/standard/upf"

    # Reset the external dependencies of the model
    model.spin_orbit = "wo_soc"

    # Test structure-dependent family change
    silicon = generate_structure_data("silicon")
    model.input_structure = silicon
    assert "Si" in model.dictionary.keys()
    assert model.ecutwfc == 30
    assert model.ecutrho == 240

    # Test that changing the structure triggers a reset
    silica = generate_structure_data("silica")
    model.input_structure = silica
    assert "Si" in model.dictionary.keys()
    assert "O" in model.dictionary.keys()

    pseudos.render()

    # Check uploaders
    assert len(pseudos.setter_widget.children) == 2

    # Check Si uploader (Si.upf)
    assert pseudos.setter_widget.children[0].kind_name == "Si"
    assert pseudos.setter_widget.children[0].kind_symbol == "Si"
    pseudo = orm.load_node(model.dictionary["Si"])
    assert pseudos.setter_widget.children[0].pseudo == pseudo
    assert pseudos.setter_widget.children[0].pseudo_text.value == pseudo.filename
    assert pseudos.setter_widget.children[0].cutoffs == [30, 240]

    # Check O uploader (O.upf)
    assert pseudos.setter_widget.children[1].kind_name == "O"
    assert pseudos.setter_widget.children[1].kind_symbol == "O"
    pseudo = orm.load_node(model.dictionary["O"])
    assert pseudos.setter_widget.children[1].pseudo == pseudo
    assert pseudos.setter_widget.children[1].pseudo_text.value == pseudo.filename
    assert pseudos.setter_widget.children[1].cutoffs == [30, 240]

    # Test reset from uploaded state
    uploader: PseudoUploadWidget = pseudos.setter_widget.children[1]
    new_O_pseudo = generate_upf_data("O", "O_new.upf", perturb=True)
    uploader._on_file_upload(
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
    pseudo = model.dictionary["O"]
    assert orm.load_node(pseudo).filename == "O_new.upf"
    assert pseudos.setter_widget.children[1].pseudo_text.value == "O_new.upf"

    model.reset()
    pseudo = model.dictionary["O"]
    assert orm.load_node(pseudo).filename != "O_new.upf"
    assert pseudos.setter_widget.children[1].pseudo_text.value != "O_new.upf"


def test_pseudo_upload_widget(generate_upf_data):
    """Test the pseudo upload widget."""

    # Test that the kind can be not the element symbol
    # the widget initialize with the pseudo as input to mock how it will
    # be used in PseudoSetter when the pseudo family is set.
    old_pseudo = generate_upf_data("O", "O_old.upf")

    w = PseudoUploadWidget(kind_name="O1", kind_symbol="O")
    w.pseudo = old_pseudo
    w.cutoffs = [30, 240]
    w.render()

    message = "ψ: <b>{ecutwfc} Ry</b> | ρ: <b>{ecutrho} Ry</b>"  # noqa: RUF001

    assert w.pseudo.filename == "O_old.upf"
    assert w.kind_name == "O1"
    assert message.format(ecutwfc=30.0, ecutrho=240.0) in w.cutoff_message.value
    assert not w.message

    # Check different element is rejected
    different_element = generate_upf_data("Si", "Si.upf")
    w._on_file_upload(
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
    assert w.pseudo.filename == "O_old.upf"
    assert "does not match" in w.message

    # Check identical content is rejected in favor of existing one
    same_content_different_name = generate_upf_data("O", "O_copy.upf")
    w._on_file_upload(
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

    assert w.pseudo.filename == "O.upf"
    assert "Identical pseudo" in w.message

    # Check different content but same filename is rejected
    different_content_same_filename = generate_upf_data("O", "O.upf", perturb=True)
    w._on_file_upload(
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
    assert w.pseudo.filename == "O.upf"
    assert "rename your file" in w.message

    # Check invalid pseudo content is rejected
    w._on_file_upload(
        {
            "new": {
                "O_invalid.upf": {
                    "content": b"<UPF version='...'>...</UPF>",
                },
            },
        }
    )
    assert "not a valid UPF file" in w.message

    # Check valid pseudo is accepted
    valid = generate_upf_data("O", "O_new.upf", perturb=True)
    w._on_file_upload(
        {
            "new": {
                "O_new.upf": {
                    "content": bytes(
                        valid.get_content(),
                        encoding="utf-8",
                    ),
                },
            },
        }
    )
    assert w.pseudo.filename == "O_new.upf"
    assert "uploaded successfully" in w.message
