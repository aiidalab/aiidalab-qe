from __future__ import annotations

import enum
import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.mixins import Confirmable, HasModels
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import WarningWidget
from aiidalab_widgets_base import LoadingWidget


class State(enum.Enum):
    """Local copy of AWB's `WizardAppWidgetStep.State`"""

    BLOCKED = -2
    FAIL = -1
    INIT = 0
    CONFIGURED = 1
    READY = 2
    ACTIVE = 3
    SUCCESS = 4


class WizardStepModel(Model):
    identifier = "qe-wizard-step"

    state = tl.UseEnum(State, default_value=State.INIT)

    def __init__(self, auto_advance: bool = True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.auto_advance = auto_advance

    @property
    def is_configured(self) -> bool:
        return self.state is State.CONFIGURED

    @property
    def is_successful(self) -> bool:
        return self.state is State.SUCCESS

    def update_state(self):
        raise NotImplementedError()


WSM = t.TypeVar("WSM", bound=WizardStepModel)


class WizardStep(ipw.VBox, t.Generic[WSM]):
    def __init__(self, model: WSM, **kwargs):
        self.loading_message = LoadingWidget(f"Loading {model.identifier} step")

        super().__init__(children=[self.loading_message], **kwargs)
        self.add_class("wizard-step")

        self._model = model

        self._model.observe(
            self._on_state_change,
            "state",
        )

        self.content = ipw.VBox()

        self.rendered = False

        self._background_class = ""

    def render(self):
        if self.rendered:
            return
        self._render()
        self.rendered = True
        self._post_render()

    def _on_state_change(self, change):
        self._update_background_color(change["new"])

    def _render(self):
        raise NotImplementedError()

    def _post_render(self):
        pass

    def _update_background_color(self, state: State):
        self.remove_class(self._background_class)
        self._background_class = f"qe-app-step-{state.name.lower()}"
        self.add_class(self._background_class)


class ConfirmableWizardStepModel(
    WizardStepModel,
    Confirmable,
):
    blocker_messages = tl.Unicode("")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.confirmation_exceptions.extend(
            [
                "state",
                "locked",
                "blockers",
                "blocker_messages",
            ]
        )

    def update_blockers(self):
        if self.state is not State.INIT:
            super().update_blockers()

    def update_blocker_messages(self):
        if self.is_blocked:
            formatted = "\n".join(f"<li>{item}</li>" for item in self.blockers)
            self.blocker_messages = f"""
                <div class="alert alert-danger" style="margin-top: 8px;">
                    <b>The step is blocked due to the following reason(s):</b>
                    <ul>
                        {formatted}
                    </ul>
                </div>
            """
        else:
            self.blocker_messages = ""

    def lock(self):
        super().lock()
        self.unobserve_all("confirmed")


CWSM = t.TypeVar("CWSM", bound=ConfirmableWizardStepModel)


class ConfirmableWizardStep(WizardStep[CWSM]):
    def __init__(
        self,
        model: CWSM,
        confirm_kwargs=None,
        **kwargs,
    ):
        super().__init__(model, **kwargs)
        self._model.observe(
            self._on_confirmation_change,
            "confirmed",
        )
        self._model.observe(
            self._on_blockers_change,
            "blockers",
        )

        if confirm_kwargs is None:
            confirm_kwargs = {}

        self.confirm_button_style = confirm_kwargs.get("button_style", "success")
        self.confirm_button_icon = confirm_kwargs.get("icon", "check-circle")
        self.confirm_button_description = confirm_kwargs.get("description", "Confirm")
        self.confirm_button_tooltip = confirm_kwargs.get("tooltip", "Confirm")

    def confirm(self, _=None):
        self._model.confirm()

    def _render(self):
        self.confirm_button = ipw.Button(
            description=self.confirm_button_description,
            tooltip=self.confirm_button_tooltip,
            button_style=self.confirm_button_style,
            icon=self.confirm_button_icon,
            layout=ipw.Layout(width="auto"),
            disabled=not self._model.is_blocked,
        )
        ipw.dlink(
            (self._model, "state"),
            (self.confirm_button, "disabled"),
            lambda _: self._model.is_blocked or not self._model.is_configured,
        )
        self.confirm_button.on_click(self.confirm)

        self.blocker_messages = ipw.HTML()
        self.blocker_messages.add_class("blocker-messages")
        ipw.dlink(
            (self._model, "blocker_messages"),
            (self.blocker_messages, "value"),
        )

        self.confirm_box = ipw.VBox(
            children=[
                self.confirm_button,
                self.blocker_messages,
            ]
        )

    def _on_confirmation_change(self, _):
        self._model.update_state()

    def _on_blockers_change(self, _):
        if self._model.state is not State.INIT:
            self._model.update_blockers()
            self._model.update_blocker_messages()
            self._model.update_state()


class DependentWizardStepModel(WizardStepModel):
    previous_step_state = tl.UseEnum(State, default_value=State.INIT)

    _dependencies: list[str] = []

    @property
    def is_previous_step_successful(self) -> bool:
        return self.previous_step_state is State.SUCCESS

    @property
    def has_all_dependencies(self) -> bool:
        return all(getattr(self, dep) for dep in self._dependencies)


DWSM = t.TypeVar("DWSM", bound=DependentWizardStepModel)


class DependentWizardStep(WizardStep[DWSM]):
    _missing_message = "Required information is missing"

    def __init__(self, model: DWSM, **kwargs):
        super().__init__(model, **kwargs)

        self.missing_message = WarningWidget(message=self._missing_message)

        self._model.observe(
            self._on_previous_step_state_change,
            "previous_step_state",
        )
        self._model.observe(
            self._update_content,
            [
                "previous_step_state",
                "state",
                *self._model._dependencies,
            ],
        )

    def _post_render(self):
        super()._post_render()
        self._update_content()

    def _on_previous_step_state_change(self, _=None):
        self._model.update_state()

    def _update_content(self, _=None):
        if not self._model.is_previous_step_successful:
            self._show_missing_info_warning()
        elif not self._model.has_all_dependencies:
            self._show_loading_message()
        else:
            self._show_content()

    def _show_missing_info_warning(self):
        self.children = [self.missing_message]

    def _show_loading_message(self):
        self.children = [self.loading_message]

    def _show_content(self):
        self.children = [self.content]


class ConfirmableDependentWizardStepModel(
    DependentWizardStepModel,
    ConfirmableWizardStepModel,
):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.confirmation_exceptions.append("previous_step_state")


CDWSM = t.TypeVar("CDWSM", bound=ConfirmableDependentWizardStepModel)


class ConfirmableDependentWizardStep(
    DependentWizardStep[CDWSM],
    ConfirmableWizardStep[CDWSM],
):
    """A confirmable dependent wizard step."""


class WizardModel(Model, HasModels[WizardStepModel]):
    state = tl.Dict(None, allow_none=True)
    selected_index = tl.Int(None, allow_none=True)
    loading = tl.Bool(False)

    def load_from_state(self, state: dict):
        pass

    def auto_advance(self):
        if self.selected_index is None:
            return

        model_list = list(self._models.values())
        index = t.cast(int, self.selected_index)

        if (
            (selected_step := model_list[index]).auto_advance
            and not (index + 1 == len(model_list))
            and selected_step.is_successful
        ):
            self.selected_index += 1


class Wizard(ipw.Accordion):
    """A wizard widget that manages multiple steps."""

    def __init__(self, model: WizardModel, icons: dict[str, str], **kwargs):
        super().__init__(**kwargs)
        self.add_class("wizard")

        self._model = model
        self._icons = icons

        self._models: list[WizardStepModel] = []
        self._steps: list[WizardStep] = []
        self._titles: dict[str, str] = {}

        ipw.link(
            (self._model, "selected_index"),
            (self, "selected_index"),
        )
        self._model.observe(
            self._on_state_change,
            "state",
        )
        self._model.observe(
            self._on_step_change,
            "selected_index",
        )

        self.rendered = False

    @property
    def current_step(self) -> int:
        return sum(
            1
            for _, model in self._model.get_models()
            if model.is_configured or model.is_successful
        )

    def add_step(
        self,
        step: WizardStep,
        model: WizardStepModel,
        title: str | None = None,
    ):
        if len(self._steps) > 0:
            previous_model = self._models[-1]
            ipw.dlink(
                (previous_model, "state"),
                (model, "previous_step_state"),
            )
        self._model.add_model(model.identifier, model)
        self._models.append(model)
        self._steps.append(step)
        self._titles[model.identifier] = title or model.identifier

    def render(self):
        if self.rendered:
            return

        self.children = self._steps
        self._update_titles()

        self.rendered = True

    def _render_step(self, step_index: int):
        step: WizardStep = self.children[step_index]  # type: ignore
        step.render()

    def _on_state_change(self, change: dict):
        self._model.load_from_state(change["new"] or {})

    def _on_step_state_change(self, _=None):
        self._update_titles()

    def _on_step_change(self, change: dict):
        if (step_index := change["new"]) is not None:
            self._render_step(step_index)

    def _update_titles(self):
        for i, (identifier, title) in enumerate(self._titles.items()):
            step_model = self._model.get_model(identifier)
            icon = self._icons.get(step_model.state, "")
            step_title = f"{icon} Step {i + 1}"
            if title:
                step_title += f": {title}"
            self.set_title(i, step_title)
