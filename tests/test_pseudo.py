def test_pseudos_family_selector_widget():
    """Test the pseudos widget."""
    from aiidalab_qe.app.configuration.pseudos import PseudoFamilySelector

    w = PseudoFamilySelector()
    w.override_protocol_pseudo_family.value = True
    w.protocol_selection.value = "PseudoDojo stringent"
    w.dft_functional.value = "PBE"
    assert w.value == "PseudoDojo/0.4/PBE/SR/stringent/upf"


def test_pseudos_setter_widget(structure_data_object, generate_upf_data):
    """Test the pseudo setter widget."""
    from aiidalab_qe.app.configuration.pseudos import PseudoSetter

    # test the widget is set with the elements of the structure
    w = PseudoSetter()
    w.structure = structure_data_object("BaTiO3")
    w.override_pseudos.value = True

    assert list(w.pseudo_setter_dict.keys()) == ["Ba", "Ti", "O"]
    assert w.pseudo_setter_dict["Ba"].value.filename == "Ba.upf"
    assert w.pseudo_setter_dict["Ti"].value.filename == "Ti.upf"

    # create a custom pseudo for Ba and set to the widget
    custum_filename = "Ba_custom.upf"
    pseudo = generate_upf_data("Ba", custum_filename)
    w.pseudo_setter_dict["Ba"].value = pseudo
    assert w.pseudo_setter_dict["Ba"].value.filname == custum_filename


def test_pseudo_upload_widget(generate_upf_data):
    """Test the pseudo upload widget."""
    from aiidalab_qe.app.pseudos import PseudoUploadWidget

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
    w.ecutrho_setter.value = 250
    w.ecutwfc_setter.value = 35

    assert w.pseudo.filename == "O_new.upf"
    assert w.kind == "O1"
    assert w.ecutrho == 250
    assert w.ecutwfc == 35
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
