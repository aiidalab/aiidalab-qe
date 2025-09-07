import typing as t

from aiidalab_qe.app.submission import SubmissionStep, SubmissionStepModel
from aiidalab_qe.app.submission.global_settings import GlobalResourceSettingsModel
from aiidalab_qe.app.wizard import Wizard
from aiidalab_qe.common.code import PwCodeModel
from aiidalab_qe.common.widgets import PwCodeResourceSetupWidget
from aiidalab_qe.common.wizard import State
from aiidalab_qe.utils import shallow_copy_nested_dict


def test_code_not_selected(submit_app_generator):
    """Test if there is an error when the code is not selected."""
    app: Wizard = submit_app_generator(properties=["dos"])
    model = app.submit_model
    model.get_model("global").get_model("quantumespresso__dos").selected = None
    # Check builder construction passes without an error
    parameters = shallow_copy_nested_dict(app.submit_model.input_parameters)
    parameters |= {"codes": app.submit_model.get_model_state()}
    model._create_builder(parameters)


def test_set_codes(submit_app_generator):
    """Test setting codes (in practice, from a loaded process)."""
    app: Wizard = submit_app_generator()
    resources = app.submit_model.get_model_state()
    model = SubmissionStepModel()
    _ = SubmissionStep(model=model, auto_setup=False)
    model.await_resources()
    for identifier, code_model in app.submit_model.get_model("global").get_models():
        model.get_model("global").get_model(identifier).is_active = code_model.is_active
    model.get_model("global").set_model_state(resources["global"])  # type: ignore
    model.previous_step_state = State.SUCCESS
    assert model.get_model_state() == app.submit_model.get_model_state()


def test_global_code_toggle(app: Wizard):
    """Test that global codes toggle on/off based on their activity."""
    app.submit_model.await_resources()

    global_resources_model = t.cast(
        GlobalResourceSettingsModel,
        app.submit_model.get_model("global"),
    )
    global_resources = app.submit_step.global_resources
    global_resources.render()

    dos_code_model = global_resources_model.get_model("quantumespresso__dos")
    assert dos_code_model.is_active is False

    global_resources_model.input_parameters = {"workchain": {"properties": ["pdos"]}}
    assert dos_code_model.is_active is True
    assert global_resources.code_widgets["dos"].layout.display == "block"

    global_resources_model.input_parameters = {"workchain": {"properties": []}}
    assert dos_code_model.is_active is False
    assert global_resources.code_widgets["dos"].layout.display == "none"


def test_check_blockers(app_to_submit: Wizard):
    """Test check_submission_blockers method."""
    model = app_to_submit.submit_model
    model.await_resources()

    assert len(model.blockers) == 0

    # set dos code to None, will introduce another blocker
    dos_code = model.get_model("global").get_model("quantumespresso__dos")
    dos_value = dos_code.selected
    dos_code.selected = None
    assert len(model.blockers) == 1

    # set dos code back will remove the blocker
    dos_code.selected = dos_value
    assert len(model.blockers) == 0

    model.input_parameters = {}
    assert len(model.blockers) == 1
    assert "input parameters" in model.blockers[0]  # type: ignore

    model.structure_uuid = None
    assert len(model.blockers) == 2
    assert "input structure" in model.blockers[0]  # type: ignore
    assert "input parameters" in model.blockers[1]  # type: ignore


def test_qeapp_computational_resources_widget(app: Wizard):
    """Test QEAppComputationalResourcesWidget."""
    app.submit_model.await_resources()
    app.submit_step.render()
    global_model = app.submit_model.get_model("global")
    global_resources = app.submit_step.global_resources
    pw_code_model = t.cast(
        PwCodeModel,
        global_model.get_model("quantumespresso__pw"),
    )
    pw_code_widget = t.cast(
        PwCodeResourceSetupWidget,
        global_resources.code_widgets["pw"],
    )
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
