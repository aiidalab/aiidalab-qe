import pytest
from aiida import orm

from aiidalab_qe.common.setup_pseudos import (
    PSEUDODOJO_VERSION,
    SSSP_VERSION,
    _construct_cmd,
    _install_pseudos,
    pseudos_to_install,
)


def test_setup_pseudos_cmd():
    """Test _construct_cmd function in setup_pseudos.py."""

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


@pytest.mark.usefixtures("aiida_profile_clean")
def test_pseudos_installation():
    """Test install_pseudos"""
    # Test by compare the pseudos_to_install before and after the installation
    assert len(pseudos_to_install()) == 8
    EXPECTED_PSEUDOS = {
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/standard/upf",
        f"SSSP/{SSSP_VERSION}/PBE/efficiency",
    }

    # Install the pseudos
    [_ for _ in _install_pseudos(EXPECTED_PSEUDOS)]

    # Two pseudos are installed
    assert len(pseudos_to_install()) == 6


@pytest.mark.usefixtures("aiida_profile_clean")
def test_download_and_install_pseudo_from_file(tmp_path):
    """Test download and install pseudo from file."""
    assert len(pseudos_to_install()) == 8
    EXPECTED_PSEUDOS = {
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/standard/upf",
        f"SSSP/{SSSP_VERSION}/PBE/efficiency",
    }

    # raise error if the file does not exist
    with pytest.raises(FileNotFoundError):
        [_ for _ in _install_pseudos(EXPECTED_PSEUDOS, cwd=tmp_path)]

    # Download the pseudos to the tmp_path but not install
    [_ for _ in _install_pseudos(EXPECTED_PSEUDOS, download_only=True, cwd=tmp_path)]

    assert len(pseudos_to_install()) == 8
    assert len(list(tmp_path.iterdir())) == 2

    # Install the pseudos from the tmp_path
    [_ for _ in _install_pseudos(EXPECTED_PSEUDOS, cwd=tmp_path)]

    # Two pseudos are installed
    assert len(pseudos_to_install()) == 6


def test_pseudos_family_selector_widget():
    """Test the pseudos widget."""
    from aiidalab_qe.app.configuration.pseudos import PseudoFamilySelector

    w = PseudoFamilySelector()
    assert w.override.value is False

    w.override.value = True

    # test the default value
    assert w.value == "SSSP/1.2/PBEsol/efficiency"

    # Test if the protocol change the value will be updated
    w.protocol = "precise"
    assert w.value == "SSSP/1.2/PBEsol/precision"

    # test the functional change will update the value
    w.dft_functional.value = "PBE"
    assert w.value == "SSSP/1.2/PBE/precision"

    # Test if selecet new pseudo library the value will be updated
    w.library_selection.value = "PseudoDojo stringent"
    assert w.value == "PseudoDojo/0.4/PBE/SR/stringent/upf"


@pytest.mark.usefixtures("sssp")
def test_pseudos_setter_widget(generate_structure_data, generate_upf_data):
    """Test the pseudo setter widget."""
    from aiidalab_qe.app.configuration.pseudos import PseudoSetter

    # test the widget is set with the elements of the structure
    silicon = generate_structure_data("silicon")
    w = PseudoSetter(structure=silicon, pseudo_family="SSSP/1.2/PBEsol/efficiency")

    assert "Si" in w.pseudos.keys()
    assert w.ecutwfc == 30
    assert w.ecutrho == 240

    # reset the structure, the widget should be reset
    silica = generate_structure_data("silica")
    w.structure = silica
    assert "Si" in w.pseudos.keys()
    assert "O" in w.pseudos.keys()

    # Upload and set a new pseudo for O
    new_O_pseudo = generate_upf_data("O", "O_new.upf")
    upload_w = w.pseudo_setting_widgets.children[1]
    upload_w._on_file_upload(
        {
            "new": {
                "O_new.upf": {
                    "content": bytes(new_O_pseudo.get_content(), encoding="utf-8"),
                },
            },
        }
    )

    assert orm.load_node(w.pseudos["O"]).filename == "O_new.upf"
    #
    pseudos = w.pseudos
    cutoffs = {"cutoff_wfc": w.ecutwfc, "cutoff_rho": w.ecutrho}
    w._reset()
    assert orm.load_node(w.pseudos["O"]).filename != "O_new.upf"
    w.set_pseudos(pseudos, cutoffs)
    assert orm.load_node(w.pseudos["O"]).filename == "O_new.upf"


def test_pseudo_upload_widget(generate_upf_data):
    """Test the pseudo upload widget."""
    from aiidalab_qe.app.configuration.pseudos import PseudoUploadWidget

    # Test that the kind can be not the element symbol
    # the widget initialize with the pseudo as input to mock how it will
    # be used in PseudoSetter when the pseudo family is set.
    old_pseudo = generate_upf_data("O", "O_old.upf")
    w = PseudoUploadWidget(
        kind="O1", pseudo=old_pseudo, cutoffs={"cutoff_wfc": 30, "cutoff_rho": 240}
    )

    # when initialized the traitlets should be set already
    assert w.pseudo.filename == "O_old.upf"
    assert w.kind == "O1"
    assert w.ecutrho == 240
    assert w.ecutwfc == 30
    assert w.error_message is None

    # mimic upload a new pseudo and set cutoffs
    new_pseudo = generate_upf_data("O", "O_new.upf")
    w._on_file_upload(
        {
            "new": {
                "O_new.upf": {
                    "content": bytes(new_pseudo.get_content(), encoding="utf-8")
                }
            }
        }
    )

    assert w.pseudo.filename == "O_new.upf"
    assert w.kind == "O1"
    assert w.error_message is None

    # test upload a invalid pseudo of other element
    invalid_pseudo = generate_upf_data("Ba", "Ba.upf")
    w._on_file_upload(
        {
            "new": {
                "Ba_upf.upf": {
                    "content": bytes(invalid_pseudo.get_content(), encoding="utf-8")
                }
            }
        }
    )

    assert w.error_message is not None
    assert w.pseudo is None
    assert w.ecutrho is None
    assert w.ecutrho is None
