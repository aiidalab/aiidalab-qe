from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.app.configuration import ConfigurationStep, ConfigurationStepModel
from aiidalab_qe.app.result import ResultsStep, ResultsStepModel
from aiidalab_qe.app.structure import StructureStep, StructureStepModel
from aiidalab_qe.app.submission import SubmissionStep, SubmissionStepModel
from aiidalab_qe.common.wizard import State, WizardStep

from .model import WizardModel


class Wizard(ipw.Accordion):
    """The main widget that combines all the application steps together."""

    ICONS = {
        State.INIT: "\u25cb",
        State.READY: "\u25ce",
        State.CONFIGURED: "\u25cf",
        State.ACTIVE: "\u231b",
        State.SUCCESS: "\u2713",
        State.FAIL: "\u00d7",
    }

    TITLES = {
        "structure": "Select structure",
        "configure": "Configure workflow",
        "submit": "Choose computational resources",
        "results": "Status & results",
    }

    def __init__(
        self,
        model: WizardModel,
        auto_setup: bool = True,
        log_widget: ipw.Output | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self._model = model
        ipw.link(
            (self._model, "selected_index"),
            (self, "selected_index"),
        )

        self.structure_model = StructureStepModel(auto_advance=True)
        self.structure_step = StructureStep(
            model=self.structure_model,
            auto_setup=auto_setup,
        )
        self._model.add_model("structure", self.structure_model)

        self.configure_model = ConfigurationStepModel(auto_advance=True)
        self.configure_step = ConfigurationStep(model=self.configure_model)
        self._model.add_model("configure", self.configure_model)
        ipw.dlink(
            (self.structure_model, "state"),
            (self.configure_model, "previous_step_state"),
        )

        self.submit_model = SubmissionStepModel(auto_advance=True)
        self.submit_step = SubmissionStep(
            model=self.submit_model,
            auto_setup=auto_setup,
        )
        self._model.add_model("submit", self.submit_model)
        ipw.dlink(
            (self.configure_model, "state"),
            (self.submit_model, "previous_step_state"),
        )

        self.results_model = ResultsStepModel()
        self.results_step = ResultsStep(
            model=self.results_model,
            log_widget=log_widget,
        )
        self._model.add_model("results", self.results_model)
        ipw.dlink(
            (self.submit_model, "state"),
            (self.results_model, "previous_step_state"),
        )

        for step_model in (
            self.structure_model,
            self.configure_model,
            self.submit_model,
            self.results_model,
        ):
            step_model.observe(
                self._on_step_state_change,
                "state",
            )

        self.structure_model.observe(
            self._on_structure_confirmation_change,
            "confirmed",
        )
        self.configure_model.observe(
            self._on_configuration_confirmation_change,
            "confirmed",
        )
        self.submit_model.observe(
            self._on_submission,
            "confirmed",
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

        self.render()

    @property
    def current_step(self) -> int:
        return sum(
            1
            for _, model in self._model.get_models()
            if model.is_configured or model.is_finished
        )

    def render(self):
        if self.rendered:
            return

        self.children = [
            self.structure_step,
            self.configure_step,
            self.submit_step,
            self.results_step,
        ]
        self._update_titles()

        self.rendered = True

    def _on_state_change(self, change: dict):
        self._model.load_from_state(change["new"] or {})

    def _on_step_state_change(self, _=None):
        self._update_titles()

    def _on_step_change(self, change: dict):
        if (step_index := change["new"]) is not None:
            self._render_step(step_index)

    def _on_structure_confirmation_change(self, _):
        self._model.auto_advance()
        self._model.update_configuration_model()

    def _on_configuration_confirmation_change(self, _):
        self._model.auto_advance()
        self._model.update_submission_model()

    def _on_submission(self, _):
        self._model.auto_advance()
        self._model.update_results_model()
        self._model.lock_app()  # TODO .lock() might be enough - check!

    def _render_step(self, step_index: int):
        step: WizardStep = self.children[step_index]  # type: ignore
        step.render()

    def _update_titles(self):
        for i, (step_id, title) in enumerate(self.TITLES.items()):
            step_model = self._model.get_model(step_id)
            icon = self.ICONS.get(step_model.state)
            self.set_title(i, f"{icon} Step {i + 1}: {title}")
