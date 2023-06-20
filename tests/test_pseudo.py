def test_pseudos_family_selector_widget():
    """Test the pseudos widget."""
    from aiidalab_qe.app.pseudos import PseudoFamilySelector

    wg = PseudoFamilySelector()
    wg.override_protocol_pseudo_family.value = True
    wg.protocol_selection.value = "PseudoDojo stringent"
    wg.dft_functional.value = "PBE"
    assert wg.value == "PseudoDojo/0.4/PBE/SR/stringent/upf"
