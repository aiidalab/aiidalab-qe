from aiida import orm


def test_pseudos_family_selector_widget_protocol():
    """Test the pseudos widget."""
    from aiidalab_qe.app.configuration.pseudos import PseudoFamilySelector

    w = PseudoFamilySelector()
    w.protocol = "precise"
    assert w.protocol_selection.value == "SSSP precision"
    assert w.value == "SSSP/1.2/PBEsol/precision"


def test_pseudos_family_selector_widget():
    """Test the pseudos widget."""
    from aiidalab_qe.app.configuration.pseudos import PseudoFamilySelector

    w = PseudoFamilySelector()
    w.override_protocol_pseudo_family.value = True
    w.protocol_selection.value = "PseudoDojo stringent"
    w.dft_functional.value = "PBE"
    assert w.value == "PseudoDojo/0.4/PBE/SR/stringent/upf"


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
