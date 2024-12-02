from aiidalab_qe.app.main import App
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_qe.app.submission.model import SubmissionStepModel


def test_code_not_selected(submit_app_generator):
    """Test if there is an error when the code is not selected."""
    app: App = submit_app_generator(properties=["dos"])
    model = app.submit_model
    model.get_model("global").get_code("quantumespresso.dos").selected = None
    # Check builder construction passes without an error
    parameters = model.get_model_state()
    model._create_builder(parameters)


def test_set_selected_codes(submit_app_generator):
    """Test set_selected_codes method."""
    app: App = submit_app_generator()
    parameters = app.submit_model.get_model_state()
    model = SubmissionStepModel()
    _ = SubmitQeAppWorkChainStep(model=model, qe_auto_setup=False)
    for name, code_model in app.submit_model.get_model("global").codes.items():
        model.get_model("global").get_code(name).is_active = code_model.is_active
    model.qe_installed = True
    model.get_model("global").set_selected_codes(parameters["codes"]["global"]["codes"])
    assert model.get_selected_codes() == app.submit_model.get_selected_codes()


def test_update_codes_display(app: App):
    """Test update_codes_display method.
    If the workchain property is not selected, the related code should be hidden.
    """
    app.submit_step.render()
    model = app.submit_model
    model.get_model("global").update_active_codes()
    assert (
        app.submit_step.global_code_settings.code_widgets["dos"].layout.display
        == "none"
    )
    model.input_parameters = {"workchain": {"properties": ["pdos"]}}
    model.get_model("global").update_active_codes()
    assert (
        app.submit_step._model.get_model("global")
        .codes["quantumespresso.dos"]
        .is_active
        is True
    )
    assert (
        app.submit_step.global_code_settings.code_widgets["dos"].layout.display
        == "block"
    )


def test_check_submission_blockers(app: App):
    """Test check_submission_blockers method."""
    model = app.submit_model

    model.update_submission_blockers()
    assert len(model.internal_submission_blockers) == 0

    model.input_parameters = {"workchain": {"properties": ["pdos"]}}
    model.update_submission_blockers()
    assert len(model.internal_submission_blockers) == 0

    # set dos code to None, will introduce another blocker
    dos_code = model.get_model("global").get_code("quantumespresso.dos")
    dos_value = dos_code.selected
    dos_code.selected = None
    model.update_submission_blockers()
    assert len(model.internal_submission_blockers) == 1

    # set dos code back will remove the blocker
    dos_code.selected = dos_value
    model.update_submission_blockers()
    assert len(model.internal_submission_blockers) == 0


def test_qeapp_computational_resources_widget(app: App):
    """Test QEAppComputationalResourcesWidget."""
    app.submit_step.render()
    pw_code_model = app.submit_model.get_model("global").get_code("quantumespresso.pw")
    pw_code_widget = app.submit_step.global_code_settings.code_widgets["pw"]
    assert pw_code_widget.parallelization.npool.layout.display == "none"
    pw_code_model.override = True
    pw_code_model.npool = 2
    assert pw_code_widget.parallelization.npool.layout.display == "block"
    assert pw_code_widget.parameters == {
        "code": app.submit_step.global_code_settings.code_widgets[
            "pw"
        ].value,  # TODO why None?
        "cpus": 1,
        "cpus_per_task": 1,
        "max_wallclock_seconds": 43200,
        "nodes": 1,
        "ntasks_per_node": 1,
        "parallelization": {"npool": 2},
    }
