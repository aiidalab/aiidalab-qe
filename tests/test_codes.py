import pytest

from aiidalab_qe.app.main import App
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_qe.app.submission.model import SubmissionModel


@pytest.mark.usefixtures("aiida_profile_clean", "sssp")
def test_code_not_selected(submit_app_generator):
    """Test if there is an error when the code is not selected."""
    app: App = submit_app_generator(properties=["dos"])
    model = app.submit_model
    model.code_widgets["dos"].value = None
    # Check builder construction passes without an error
    parameters = model._get_submission_parameters()
    model._create_builder(parameters)


@pytest.mark.usefixtures("aiida_profile_clean", "sssp")
def test_set_selected_codes(submit_app_generator):
    """Test set_selected_codes method."""
    app: App = submit_app_generator()
    parameters = app.submit_model._get_submission_parameters()
    model = SubmissionModel()
    new_submit_step = SubmitQeAppWorkChainStep(model=model, qe_auto_setup=False)
    new_submit_step.render()
    model.set_selected_codes(parameters["codes"])
    assert model.get_selected_codes() == app.submit_model.get_selected_codes()


def test_update_codes_display(app: App):
    """Test update_codes_display method.
    If the workchain property is not selected, the related code should be hidden.
    """
    model = app.submit_model
    model.update_active_codes()
    assert model.code_widgets["dos"].layout.display == "none"
    model.input_parameters = {"workchain": {"properties": ["pdos"]}}
    model.update_active_codes()
    assert model.code_widgets["dos"].layout.display == "block"


@pytest.mark.usefixtures("aiida_profile_clean", "sssp")
def test_check_submission_blockers(app: App):
    """Test check_submission_blockers method."""
    model = app.submit_model

    blockers = list(model._check_submission_blockers())
    assert len(blockers) == 0

    model.input_parameters = {"workchain": {"properties": ["pdos"]}}
    blockers = list(model._check_submission_blockers())
    assert len(blockers) == 0

    # set dos code to None, will introduce another blocker
    dos_code_widget = model.code_widgets["dos"]
    dos_value = dos_code_widget.value
    dos_code_widget.value = None
    blockers = list(model._check_submission_blockers())
    assert len(blockers) == 1

    # set dos code back will remove the blocker
    dos_code_widget.value = dos_value
    blockers = list(model._check_submission_blockers())
    assert len(blockers) == 0


def test_qeapp_computational_resources_widget(app: App):
    """Test QEAppComputationalResourcesWidget."""

    pw_code_widget = app.submit_model.code_widgets["pw"]
    assert pw_code_widget.parallelization.npool.layout.display == "none"
    pw_code_widget.parallelization.override.value = True
    pw_code_widget.parallelization.npool.value = 2
    assert pw_code_widget.parallelization.npool.layout.display == "block"
    assert pw_code_widget.parameters == {
        "code": app.submit_model.code_widgets["pw"].value,  # TODO why None?
        "cpus": 1,
        "cpus_per_task": 1,
        "max_wallclock_seconds": 43200,
        "nodes": 1,
        "ntasks_per_node": 1,
        "parallelization": {"npool": 2},
    }
