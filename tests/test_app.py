import pytest

from aiidalab_qe.app.main import App


@pytest.mark.usefixtures("aiida_profile_clean", "sssp")
def test_reload_and_reset(generate_qeapp_workchain):
    app = App(qe_auto_setup=False)
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
    from aiidalab_qe.app.structure.model import StructureModel

    model = StructureModel()
    model.structure = generate_structure_data()
    assert model.structure is not None
    model.confirm()
    model.structure = generate_structure_data()
    assert not model.confirmed


@pytest.mark.usefixtures("aiida_profile_clean", "sssp")
def test_unsaved_changes(app_to_submit):
    """Test if the unsaved changes are handled correctly"""
    from aiidalab_widgets_base import WizardAppWidgetStep

    app: App = app_to_submit
    # go to the configure step, and make some changes
    app._wizard_app_widget.selected_index = 1
    app.configure_model.relax_type = "positions"
    # go to the submit step
    app._wizard_app_widget.selected_index = 2
    # the state of the configure step should be updated.
    assert app.configure_step.state == WizardAppWidgetStep.State.CONFIGURED
    # check if a new blocker is added
    assert len(app.submit_model.external_submission_blockers) == 1
    # confirm the changes
    app._wizard_app_widget.selected_index = 1
    app.configure_step.confirm()
    app._wizard_app_widget.selected_index = 2
    # the blocker should be removed
    assert len(app.submit_model.external_submission_blockers) == 0
