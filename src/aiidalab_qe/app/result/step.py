from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.wizard import DependentWizardStep, State
from aiidalab_widgets_base import LoadingWidget, ProcessMonitor

from .components import ResultsComponent
from .components.status import WorkChainStatusModel, WorkChainStatusPanel
from .components.summary import WorkflowSummary, WorkflowSummaryModel
from .components.viewer import ResultsViewer, ResultsViewerModel
from .model import ResultsStepModel


class ResultsStep(DependentWizardStep[ResultsStepModel]):
    def __init__(
        self,
        model: ResultsStepModel,
        log_widget: ipw.Output | None,
        **kwargs,
    ):
        super().__init__(model=model, **kwargs)

        self.log_widget = log_widget

        summary_model = WorkflowSummaryModel()
        self.summary_panel = WorkflowSummary(model=summary_model)
        self._model.add_model("summary", summary_model)

        results_model = ResultsViewerModel()
        self.results_panel = ResultsViewer(model=results_model)
        self._model.add_model("results", results_model)

        status_model = WorkChainStatusModel()
        self.status_panel = WorkChainStatusPanel(model=status_model)
        self._model.add_model("status", status_model)

        self.panels = {
            "Summary": self.summary_panel,
            "Status": self.status_panel,
            "Results": self.results_panel,
        }

        self._model.observe(
            self._on_previous_step_state_change,
            "previous_step_state",
        )
        self._model.observe(
            self._on_process_change,
            "process_uuid",
        )

    def reset(self):
        self._model.reset()

    def _render(self):
        self.kill_button = ipw.Button(
            description="Kill workchain",
            tooltip="Kill the below workchain.",
            button_style="danger",
            icon="stop",
            layout=ipw.Layout(width="auto", display="none"),
        )
        ipw.dlink(
            (self._model, "state"),
            (self.kill_button, "disabled"),
            lambda state: state is not State.ACTIVE,
        )
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.kill_button.layout, "display"),
            lambda _: "none"
            if not self._model.has_process
            or self._model.process.is_finished
            or self._model.process.is_excepted
            or self._model.state in {State.SUCCESS, State.FAIL}
            else "block",
        )
        self.kill_button.on_click(self._on_kill_button_click)

        self.clean_scratch_button = ipw.Button(
            description="Clean remote data",
            tooltip="Clean the remote folders of the workchain.",
            button_style="danger",
            icon="trash",
            layout=ipw.Layout(width="auto", display="none"),
        )
        ipw.dlink(
            (self._model, "process_remote_folder_is_clean"),
            (self.clean_scratch_button, "disabled"),
        )
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.clean_scratch_button.layout, "display"),
            lambda _: "block"
            if self._model.has_process and self._model.process.is_terminated
            else "none",
        )
        self.clean_scratch_button.on_click(self._on_clean_scratch_button_click)

        self.process_info = ipw.HTML()
        ipw.dlink(
            (self._model, "process_info"),
            (self.process_info, "value"),
        )

        self.toggle_controls = ipw.ToggleButtons(
            options=[*self.panels.keys()],
            tooltips=[
                "A summary of calculation parameters",
                "The calculation results",
                "A detailed progress status of the workflow",
            ],
            icons=[
                "file-text-o",
                "bar-chart",
                "tasks",
            ],
            value=None,
        )
        self.toggle_controls.add_class("results-step-toggles")
        self.toggle_controls.observe(
            self._on_toggle_change,
            "value",
        )

        self.container = ipw.VBox(
            children=[
                self.status_panel,
            ],
        )

        loading_message = LoadingWidget(message="Loading results")

        ipw.dlink(
            (self._model, "process_uuid"),
            (self, "children"),
            lambda _: (
                [
                    InAppGuide(identifier="results-step"),
                    self.process_info,
                    ipw.HBox(
                        children=[
                            self.kill_button,
                            self.clean_scratch_button,
                        ],
                        layout=ipw.Layout(margin="0 3px"),
                    ),
                    self.toggle_controls,
                    self.container,
                ]
                if self._model.has_process
                else [loading_message]
            )
            if self._model.is_previous_step_successful
            else [self._model.missing_process_warning],
        )

    def _post_render(self):
        self.toggle_controls.value = (
            "Results"
            if self._model.has_process and self._model.process.is_finished_ok
            else "Status"
        )

        self.process_monitor = ProcessMonitor(
            timeout=0.5,
            callbacks=[
                self._update_status,
                lambda _: self._model.update_state(),
            ],
            log_widget=self.log_widget,
        )
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.process_monitor, "value"),
        )

    def _on_previous_step_state_change(self, _):
        if self._model.is_previous_step_successful:
            message = (
                "Loading results"
                if self._model.has_process and self._model.process.is_finished
                else "Submitting calculation"
            )
            self.children = [LoadingWidget(message)]

    def _on_toggle_change(self, change):
        panel = self.panels[change["new"]]
        self._toggle_view(panel)

    def _on_process_change(self, _):
        self._model.update()
        self._model.update_state()

    def _on_kill_button_click(self, _):
        self._model.kill_process()

    def _on_clean_scratch_button_click(self, _):
        self._model.clean_remote_data()

    def _toggle_view(self, panel: ResultsComponent):
        self.container.children = [panel]
        panel.render()

    def _update_status(self):
        self._model.monitor_counter += 1
