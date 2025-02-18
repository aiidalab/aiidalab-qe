from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_qe.app.submission.model import SubmissionStepModel
from aiidalab_qe.app.wizard_app import WizardApp


def test_code_not_selected(submit_app_generator):
    """Test if there is an error when the code is not selected."""
    app: WizardApp = submit_app_generator(properties=["dos"])
    model = app.submit_model
    model.get_model("global").get_model("quantumespresso__dos").selected = None
    # Check builder construction passes without an error
    parameters = model.get_model_state()
    model._create_builder(parameters)


def test_set_selected_codes(submit_app_generator):
    """Test set_selected_codes method."""
    app: WizardApp = submit_app_generator()
    parameters = app.submit_model.get_model_state()
    model = SubmissionStepModel()
    _ = SubmitQeAppWorkChainStep(model=model, auto_setup=False)
    for identifier, code_model in app.submit_model.get_model("global").get_models():
        model.get_model("global").get_model(identifier).is_active = code_model.is_active
    model.qe_installed = True
    model.get_model("global").set_selected_codes(parameters["codes"]["global"]["codes"])
    assert model.get_selected_codes() == app.submit_model.get_selected_codes()


def test_update_codes_display(app: WizardApp):
    """Test update_codes_display method.
    If the workchain property is not selected, the related code should be hidden.
    """
    app.submit_step.render()
    model = app.submit_model
    global_model = model.get_model("global")
    global_model.update_active_codes()
    global_resources = app.submit_step.global_resources
    assert global_resources.code_widgets["dos"].layout.display == "none"
    model.input_parameters = {"workchain": {"properties": ["pdos"]}}
    global_model.update_active_codes()
    assert global_model.get_model("quantumespresso__dos").is_active is True
    assert global_resources.code_widgets["dos"].layout.display == "block"


def test_check_blockers(app: WizardApp):
    """Test check_submission_blockers method."""
    model = app.submit_model

    model.update_blockers()
    assert len(model.blockers) == 0

    model.input_parameters = {"workchain": {"properties": ["pdos"]}}
    model.update_blockers()
    assert len(model.blockers) == 0

    # set dos code to None, will introduce another blocker
    dos_code = model.get_model("global").get_model("quantumespresso__dos")
    dos_value = dos_code.selected
    dos_code.selected = None
    model.update_blockers()
    assert len(model.blockers) == 1

    # set dos code back will remove the blocker
    dos_code.selected = dos_value
    model.update_blockers()
    assert len(model.blockers) == 0


def test_qeapp_computational_resources_widget(app: WizardApp):
    """Test QEAppComputationalResourcesWidget."""
    app.submit_step.render()
    global_model = app.submit_model.get_model("global")
    global_resources = app.submit_step.global_resources
    pw_code_model = global_model.get_model("quantumespresso__pw")
    pw_code_widget = global_resources.code_widgets["pw"]
    assert pw_code_widget.parallelization.npool.layout.display == "none"
    pw_code_model.parallelization_override = True
    pw_code_model.npool = 2
    assert pw_code_widget.parallelization.npool.layout.display == "block"
    assert pw_code_widget.parameters == {
        "code": global_resources.code_widgets["pw"].value,
        "cpus": 1,
        "cpus_per_task": 1,
        "max_wallclock_seconds": 43200,
        "nodes": 1,
        "ntasks_per_node": 1,
        "parallelization": {"npool": 2},
    }
