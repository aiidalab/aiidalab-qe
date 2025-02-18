from aiidalab_qe.app.wizard_app import WizardApp


def test_reload_and_reset(generate_qeapp_workchain):
    app = WizardApp(auto_setup=False)
    workchain = generate_qeapp_workchain(
        relax_type="positions",
        spin_type="collinear",
        run_bands=True,
        run_pdos=False,
    )

    # Test if the app can be loaded from process
    app.process = workchain.node.pk
    assert app.configure_model.relax_type == "positions"
    assert app.configure_model.get_model("workchain").spin_type == "collinear"
    assert app.configure_model.get_model("bands").include is True
    assert app.configure_model.get_model("pdos").include is False
    advanced_model = app.configure_model.get_model("advanced")
    assert len(advanced_model.get_model("pseudos").dictionary) > 0
    assert app.configure_step.state == app.configure_step.State.SUCCESS


def test_selecting_new_structure_unconfirms_model(generate_structure_data):
    from aiidalab_qe.app.structure.model import StructureStepModel

    model = StructureStepModel()
    model.input_structure = generate_structure_data()
    assert model.has_structure
    model.confirm()
    model.input_structure = generate_structure_data()
    assert not model.confirmed
