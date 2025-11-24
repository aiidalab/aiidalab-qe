from __future__ import annotations

import enum
import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.mixins import Confirmable, HasBlockers
from aiidalab_qe.common.mvc import Model
from aiidalab_widgets_base import LoadingWidget


class State(enum.Enum):
    """Local copy of AWB's `WizardAppWidgetStep.State`"""

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

        self._model = model

        self._model.observe(
            self._on_state_change,
            "state",
        )

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
    HasBlockers,
):
    blockers = tl.List(tl.Unicode())
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

    @property
    def is_blocked(self):
        return any(self.blockers)

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
        self.content = ipw.VBox()

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
        self.children += (self.confirm_box,)

    def _on_confirmation_change(self, _):
        self._model.update_state()

    def _on_blockers_change(self, _):
        self._model.update_blocker_messages()
        self._model.update_state()


class DependentWizardStepModel(
    WizardStepModel,
):
    previous_step_state = tl.UseEnum(State, default_value=State.INIT)

    @property
    def is_previous_step_successful(self) -> bool:
        return self.previous_step_state is State.SUCCESS


DWSM = t.TypeVar("DWSM", bound=DependentWizardStepModel)


class DependentWizardStep(WizardStep[DWSM]):
    def __init__(self, model: DWSM, **kwargs):
        super().__init__(model, **kwargs)
        self._model.observe(
            self._on_previous_step_state_change,
            "previous_step_state",
        )

    def _on_previous_step_state_change(self, _):
        self._model.update_state()


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
