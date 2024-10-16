import pytest

from aiidalab_qe.app.main import App


# TODO: (jusong.yu) I have to add this fixture after I change sssp fixture from session to function
# The fixtures are a bit messy, we need to clean it up
# Same for other tests that use sssp fixture
@pytest.mark.usefixtures("sssp")
def test_reload_and_reset(submit_app_generator, generate_qeapp_workchain):
    """Test if the GUI paramters can be reload and reset properly"""
    workchain = generate_qeapp_workchain(
        relax_type="positions",
        run_bands=True,
        run_pdos=False,
        spin_type="collinear",
    )
    app: App = submit_app_generator()
    # select the pk
    app.work_chain_selector.value = workchain.node.pk
    # check if the value are reload correctly
    assert app.config_model.workchain.relax_type == "positions"
    assert app.config_model.workchain.spin_type == "collinear"
    assert app.config_model.get_model("bands").include is True
    assert app.config_model.get_model("pdos").include is False
    assert len(app.config_model.advanced.pseudos.dictionary) > 0
    assert app.configure_step.state == app.configure_step.State.SUCCESS
    # in the reload case, go to the submit step should not
    # trigger the reset of previous steps
    app._wizard_app_widget.selected_index = 2
    assert app.configure_step.state == app.configure_step.State.SUCCESS
    # new workflow, this will reset the GUI
    app.work_chain_selector.value = None
    # check if the value are reload correctly
    assert app.structure_step.manager.structure is None
    assert app.config_model.workchain.relax_type == "positions_cell"
    assert app.config_model.workchain.spin_type == "none"
    assert app.config_model.get_model("bands").include is False
    assert app.config_model.get_model("pdos").include is False
    assert app.config_model.advanced.pseudos.ecutwfc == 30  # TODO why 0?
    assert len(app.config_model.advanced.pseudos.dictionary) == 1  # TODO why 0?
    assert app.submit_model.code_widgets["pw"].num_cpus.value == 2  # TODO why 4?


@pytest.mark.usefixtures("sssp")
def test_select_new_structure(app_to_submit, generate_structure_data):
    """Test if the new structure is selected, the confirmed structure is reset"""
    app: App = app_to_submit
    assert app.struct_model.confirmed_structure is not None
    # select a new structure will reset the confirmed structure
    app.struct_model.structure = generate_structure_data()
    assert app.struct_model.confirmed_structure is None


@pytest.mark.usefixtures("sssp")
def test_unsaved_changes(app_to_submit):
    """Test if the unsaved changes are handled correctly"""
    from aiidalab_widgets_base import WizardAppWidgetStep

    app: App = app_to_submit
    # go to the configure step, and make some changes
    app._wizard_app_widget.selected_index = 1
    app.config_model.workchain.relax_type = "positions"
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
