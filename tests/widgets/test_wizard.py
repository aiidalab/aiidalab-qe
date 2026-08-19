import ipywidgets as ipw

from aiidalab_qe.common.wizard import (
    ConfirmableWizardStepModel,
    DependentWizardStep,
    DependentWizardStepModel,
    State,
    Wizard,
    WizardModel,
    WizardStep,
    WizardStepModel,
)


class DummyStepModel(WizardStepModel):
    def __init__(self, identifier: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.identifier = identifier

    def update_state(self):
        pass


class DummyStep(WizardStep[DummyStepModel]):
    def _render(self):
        self.content = ipw.Label("dummy content")


class DummyDependentModel(DependentWizardStepModel):
    _dependencies = ["ready"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ready = False

    def update_state(self):
        if self.previous_step_state is State.SUCCESS and self.ready:
            self.state = State.SUCCESS
        else:
            self.state = State.INIT


class DummyDependentStep(DependentWizardStep[DummyDependentModel]):
    def _render(self):
        self.content = ipw.Label("dependent content")


class DummyConfirmableModel(ConfirmableWizardStepModel):
    def update_state(self):
        if self.confirmed:
            self.state = State.SUCCESS
        else:
            self.state = State.INIT


def test_wizard_auto_advance_between_successful_steps():
    model = WizardModel()

    first = DummyStepModel("first")
    second = DummyStepModel("second")
    third = DummyStepModel("third")

    model.add_model(first.identifier, first)
    model.add_model(second.identifier, second)
    model.add_model(third.identifier, third)

    model.selected_index = 0
    first.state = State.SUCCESS

    model.auto_advance()
    assert model.selected_index == 1

    model.selected_index = 1
    second.state = State.SUCCESS

    model.auto_advance()
    assert model.selected_index == 2


def test_dependent_step_shows_warning_until_dependencies_are_met():
    model = DummyDependentModel()
    step = DummyDependentStep(model)
    step.render()

    assert step.children[0] is step.missing_message

    model.previous_step_state = State.SUCCESS
    model.ready = True
    step._update_content()
    assert step.children[0] is step.content


def test_confirmable_model_keeps_confirmation_when_state_changes():
    model = DummyConfirmableModel()

    model.confirm()
    assert model.confirmed is True

    model.state = State.SUCCESS
    assert model.confirmed is True


def test_wizard_sets_titles_and_tracks_current_step_progress():
    model = WizardModel()
    wizard = Wizard(model, icons={})

    first = DummyStepModel("first")
    second = DummyDependentModel()
    second.identifier = "second"

    wizard.add_step(DummyStep(first), first, title="First")
    wizard.add_step(DummyDependentStep(second), second, title="Second")
    wizard.render()

    assert wizard.get_title(0).endswith("First")
    assert wizard.get_title(1).endswith("Second")

    first.state = State.SUCCESS
    second.previous_step_state = State.SUCCESS
    second.state = State.CONFIGURED
    assert wizard.current_step == 2
