from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.app.configuration import ConfigurationStep, ConfigurationStepModel
from aiidalab_qe.app.result import ResultsStep, ResultsStepModel
from aiidalab_qe.app.structure import StructureStep, StructureStepModel
from aiidalab_qe.app.submission import SubmissionStep, SubmissionStepModel
from aiidalab_qe.common.wizard import State, Wizard

from .model import QeWizardModel

ICONS = {
    State.INIT: "\u25cb",
    State.READY: "\u25ce",
    State.CONFIGURED: "\u25cf",
    State.ACTIVE: "\u231b",
    State.SUCCESS: "\u2713",
    State.FAIL: "\u00d7",
    State.BLOCKED: "\u25cf",
}


class QeWizard(Wizard):
    """The main widget that combines all the application steps together."""

    def __init__(
        self,
        model: QeWizardModel,
        auto_setup: bool = True,
        log_widget: ipw.Output | None = None,
        **kwargs,
    ):
        super().__init__(model, ICONS, **kwargs)

        self.structure_model = StructureStepModel(auto_advance=True)
        self.structure_step = StructureStep(
            model=self.structure_model,
            auto_setup=auto_setup,
        )
        self.add_step(
            step=self.structure_step,
            model=self.structure_model,
            title="Select structure",
        )

        self.configuration_model = ConfigurationStepModel(auto_advance=True)
        self.configuration_step = ConfigurationStep(model=self.configuration_model)
        self.add_step(
            step=self.configuration_step,
            model=self.configuration_model,
            title="Configure workflow",
        )

        self.submission_model = SubmissionStepModel(auto_advance=True)
        self.submission_step = SubmissionStep(
            model=self.submission_model,
            auto_setup=auto_setup,
        )
        self.add_step(
            step=self.submission_step,
            model=self.submission_model,
            title="Choose computational resources",
        )

        self.results_model = ResultsStepModel()
        self.results_step = ResultsStep(
            model=self.results_model,
            log_widget=log_widget,
        )
        self.add_step(
            step=self.results_step,
            model=self.results_model,
            title="Status & results",
        )

        for model in self._models:
            model.observe(
                self._on_step_state_change,
                "state",
            )

        self.structure_model.observe(
            self._on_structure_confirmation_change,
            "confirmed",
        )
        self.configuration_model.observe(
            self._on_configuration_confirmation_change,
            "confirmed",
        )
        self.submission_model.observe(
            self._on_submission,
            "confirmed",
        )

        self.render()

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
