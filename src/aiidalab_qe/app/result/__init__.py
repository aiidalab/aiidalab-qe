import ipywidgets as ipw
import traitlets as tl

from aiida.engine import ProcessState
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.process import STATE_ICONS
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_qe.common.wizard import QeDependentWizardStep
from aiidalab_widgets_base import ProcessMonitor, WizardAppWidgetStep

from .components import ResultsComponent
from .components.status import WorkChainStatusModel, WorkChainStatusPanel
from .components.summary import WorkChainSummary, WorkChainSummaryModel
from .components.viewer import WorkChainResultsViewer, WorkChainResultsViewerModel
from .model import ResultsStepModel


class ViewQeAppWorkChainStatusAndResultsStep(QeDependentWizardStep[ResultsStepModel]):
    missing_information_warning = (
        "No available results. Did you submit or load a calculation?"
    )

    STATUS_TEMPLATE = "<h4>Workflow status: {}</h4"

    def __init__(self, model: ResultsStepModel, **kwargs):
        self.log_widget = kwargs.pop("log_widget", None)
        super().__init__(model=model, **kwargs)
        self.observe(
            self._on_previous_step_state_change,
            "previous_step_state",
        )
        self._model.observe(
            self._on_process_change,
            "process_uuid",
        )

    def _render(self):
        self.kill_button = ipw.Button(
            description="Kill workchain",
            tooltip="Kill the below workchain.",
            button_style="danger",
            icon="stop",
            layout=ipw.Layout(width="auto", display="none"),
        )
        ipw.dlink(
            (self, "state"),
            (self.kill_button, "disabled"),
            lambda state: state is not self.State.ACTIVE,
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
        self.clean_scratch_button.on_click(self._on_clean_scratch_button_click)

        self.process_info = ipw.HTML()
        ipw.dlink(
            (self._model, "process_info"),
            (self.process_info, "value"),
        )

        summary_model = WorkChainSummaryModel()
        self.summary_panel = WorkChainSummary(model=summary_model)
        self._model.add_model("summary", summary_model)

        results_model = WorkChainResultsViewerModel()
        self.results_panel = WorkChainResultsViewer(model=results_model)
        self._model.add_model("results", results_model)

        status_model = WorkChainStatusModel()
        self.status_panel = WorkChainStatusPanel(model=status_model)
        self._model.add_model("status", status_model)

        self.panels = {
            "Summary": self.summary_panel,
            "Status": self.status_panel,
            "Results": self.results_panel,
        }

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

        if self._model.has_process:
            self._update_children()

    def _post_render(self):
        self._update_kill_button_layout()
        self._update_clean_scratch_button_layout()

        self.toggle_controls.value = (
            "Results"
            if (process := self._model.fetch_process_node()) and process.is_finished_ok
            else "Status"
        )

        self.process_monitor = ProcessMonitor(
            timeout=0.5,
            callbacks=[
                self._update_status,
                self._update_state,
            ],
            log_widget=self.log_widget,
        )
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.process_monitor, "value"),
        )

    def can_reset(self):
        "Checks if process is running (active), which disallows a reset."
        return self.state is not self.State.ACTIVE

    def reset(self):
        self._model.reset()

    @tl.observe("state")
    def _on_state_change(self, change):
        super()._on_state_change(change)
        self._update_kill_button_layout()

    def _on_previous_step_state_change(self, _):
        if self.previous_step_state is WizardAppWidgetStep.State.SUCCESS:
            process_node = self._model.fetch_process_node()
            message = (
                "Loading results"
                if process_node and process_node.is_finished
                else "Submitting calculation"
            )
            self.children = [LoadingWidget(message)]

    def _on_toggle_change(self, change):
        panel = self.panels[change["new"]]
        self._toggle_view(panel)

    def _on_process_change(self, _):
        if self.rendered:
            self._update_children()
        self._model.update()
        self._update_state()
        self._update_kill_button_layout()
        self._update_clean_scratch_button_layout()

    def _on_kill_button_click(self, _):
        self._model.kill_process()
        self._update_kill_button_layout()

    def _on_clean_scratch_button_click(self, _):
        self._model.clean_remote_data()
        self._update_clean_scratch_button_layout()

    def _update_children(self):
        self.children = [
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
        self.previous_children = list(self.children)

    def _toggle_view(self, panel: ResultsComponent):
        self.container.children = [panel]
        panel.render()

    def _update_kill_button_layout(self):
        if not self.rendered:
            return
        process_node = self._model.fetch_process_node()
        if (
            not process_node
            or process_node.is_finished
            or process_node.is_excepted
            or self.state
            in (
                self.State.SUCCESS,
                self.State.FAIL,
            )
        ):
            self.kill_button.layout.display = "none"
        else:
            self.kill_button.layout.display = "block"

    def _update_clean_scratch_button_layout(self):
        if not self.rendered:
            return
        process_node = self._model.fetch_process_node()
        if process_node and process_node.is_terminated:
            self.clean_scratch_button.layout.display = "block"
        else:
            self.clean_scratch_button.layout.display = "none"

    def _update_status(self):
        self._model.monitor_counter += 1

    def _update_state(self):
        if not (process_node := self._model.fetch_process_node()):
            self.state = self.State.INIT
            self._update_controls()
            return

        if process_state := process_node.process_state:
            status = self._get_process_status(process_state.value)
        else:
            status = "Unknown"

        if process_state is ProcessState.CREATED:
            self.state = self.State.ACTIVE
        elif process_state in (
            ProcessState.RUNNING,
            ProcessState.WAITING,
        ):
            self.state = self.State.ACTIVE
            status = self._get_process_status("running")  # overwrite status
        elif process_state in (
            ProcessState.EXCEPTED,
            ProcessState.KILLED,
        ):
            self.state = self.State.FAIL
        elif process_node.is_failed:
            self.state = self.State.FAIL
        elif process_node.is_finished_ok:
            self.state = self.State.SUCCESS

        self._model.process_info = self.STATUS_TEMPLATE.format(status)

        self._update_controls()

    def _update_controls(self):
        if self.state in (self.State.SUCCESS, self.State.FAIL):
            self._update_kill_button_layout()
            self._update_clean_scratch_button_layout()

    def _get_process_status(self, state: str):
        return f"{state.capitalize()} {STATE_ICONS[state]}"
