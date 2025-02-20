from __future__ import annotations

import os
import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.mixins import Confirmable, HasBlockers, HasModels
from aiidalab_qe.common.mvc import Model
from aiidalab_widgets_base import WizardAppWidgetStep

from .widgets import LoadingWidget


class QeWizardStepModel(Model):
    identifier = "QE wizard"


WSM = t.TypeVar("WSM", bound=QeWizardStepModel)


class QeWizardStep(ipw.VBox, WizardAppWidgetStep, t.Generic[WSM]):
    def __init__(self, model: WSM, **kwargs):
        self.loading_message = LoadingWidget(f"Loading {model.identifier} step")
        super().__init__(children=[self.loading_message], **kwargs)
        self._model = model
        self.rendered = False
        self._background_class = ""

    def render(self):
        if self.rendered:
            return
        self._render()
        self.rendered = True
        self._post_render()

    @tl.observe("state")
    def _on_state_change(self, change):
        self._update_background_color(change["new"])

    def _render(self):
        raise NotImplementedError()

    def _post_render(self):
        pass

    def _update_background_color(self, state: WizardAppWidgetStep.State):
        self.remove_class(self._background_class)
        self._background_class = f"qe-app-step-{state.name.lower()}"
        self.add_class(self._background_class)

    def _update_state(self):
        pass


class QeConfirmableWizardStepModel(
    QeWizardStepModel,
    Confirmable,
    HasBlockers,
):
    blockers = tl.List(tl.Unicode())
    blocker_messages = tl.Unicode("")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.confirmation_exceptions += [
            "blockers",
            "blocker_messages",
        ]

    @property
    def is_blocked(self):
        return any(self.blockers)

    def update_blockers(self):
        blockers = list(self._check_blockers())
        if isinstance(self, HasModels):
            for _, model in self.get_models():
                if isinstance(model, HasBlockers):
                    blockers += model.blockers
        self.blockers = blockers

    def update_blocker_messages(self):
        if self.is_blocked:
            formatted = "\n".join(f"<li>{item}</li>" for item in self.blockers)
            self.blocker_messages = f"""
                <div class="alert alert-danger">
                    <b>The step is blocked due to the following reason(s):</b>
                    <ul>
                        {formatted}
                    </ul>
                </div>
            """
        else:
            self.blocker_messages = ""

    def _check_blockers(self):
        raise NotImplementedError


CWSM = t.TypeVar("CWSM", bound=QeConfirmableWizardStepModel)


class QeConfirmableWizardStep(QeWizardStep[CWSM]):
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
            (self, "state"),
            (self.confirm_button, "disabled"),
            lambda state: self._model.is_blocked or state != self.State.CONFIGURED,
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
        self._update_state()

    def _on_blockers_change(self, _):
        if self.rendered:
            self._enable_confirm_button()
        self._model.update_blocker_messages()
        self._update_state()

    def confirm(self, _=None):
        self._model.confirm()

    def _enable_confirm_button(self):
        can_confirm = self._model.is_blocked or self.state != self.State.CONFIGURED
        self.confirm_button.disabled = can_confirm


class QeDependentWizardStep(QeWizardStep[WSM]):
    missing_information_warning = "Missing information"

    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    def __init__(self, model: WSM, **kwargs):
        super().__init__(model, **kwargs)
        self.previous_children = list(self.children)
        self.warning_message = ipw.HTML(
            f"""
            <div class="alert alert-danger">
                <b>Warning:</b> {self.missing_information_warning}
            </div>
        """
        )

    def render(self):
        if "PYTEST_CURRENT_TEST" in os.environ:
            super().render()
            return
        if self.previous_step_state is WizardAppWidgetStep.State.SUCCESS:
            self._hide_missing_information_warning()
            if not self.rendered:
                super().render()
                self.previous_children = list(self.children)
        else:
            self._show_missing_information_warning()

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update_state()

    def _show_missing_information_warning(self):
        self.children = [self.warning_message]
        self.rendered = False

    def _hide_missing_information_warning(self):
        self.children = self.previous_children


class QeConfirmableDependentWizardStep(
    QeDependentWizardStep[CWSM],
    QeConfirmableWizardStep[CWSM],
):
    pass
