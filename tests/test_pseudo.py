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
