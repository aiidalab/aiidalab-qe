import ipywidgets as ipw
import traitlets as tl

from aiida.engine import ProcessState
from aiidalab_widgets_base import (
    AiidaNodeViewWidget,
    ProcessMonitor,
    ProcessNodesTreeWidget,
    WizardAppWidgetStep,
)

from .model import ResultsModel

# trigger registration of the viewer widget:
from .workchain_viewer import WorkChainViewer  # noqa: F401

PROCESS_COMPLETED = "<h4>Workflow completed successfully!</h4>"
PROCESS_EXCEPTED = "<h4>Workflow is excepted!</h4>"
PROCESS_RUNNING = "<h4>Workflow is running!</h4>"


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppWidgetStep):
    def __init__(self, model: ResultsModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading results panel")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_process_change,
            "process",
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.process_tree = ProcessNodesTreeWidget()
        ipw.dlink(
            (self._model, "process"),
            (self.process_tree, "value"),
        )

        self.node_view = AiidaNodeViewWidget(layout={"width": "auto", "height": "auto"})
        ipw.dlink(
            (self.process_tree, "selected_nodes"),
            (self.node_view, "node"),
            transform=lambda nodes: nodes[0] if nodes else None,
        )

        self.process_status = ipw.VBox(
            children=[
                self.process_tree,
                self.node_view,
            ],
        )

        self.process_monitor = ProcessMonitor(
            timeout=0.2,
            callbacks=[
                self.process_tree.update,
                self._update_state,
            ],
        )
        ipw.dlink(
            (self._model, "process"),
            (self.process_monitor, "value"),
        )

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

        self.update_results_button = ipw.Button(
            description="Update results",
            tooltip="Trigger the update of the results.",
            button_style="success",
            icon="refresh",
            layout=ipw.Layout(width="auto", display="block"),
        )
        self.update_results_button.on_click(self._on_update_results_button_click)

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

        self.children = [
            self.process_info,
            ipw.HBox(
                children=[
                    self.kill_button,
                    self.update_results_button,
                    self.clean_scratch_button,
                ],
                layout=ipw.Layout(margin="0 3px"),
            ),
            self.process_status,
        ]

        self.rendered = True

        self._update_kill_button_layout()
        self._update_clean_scratch_button_layout()

    def can_reset(self):
        return self.state is not self.State.ACTIVE

    def reset(self):
        self._model.reset()

    @tl.observe("state")
    def _on_state_change(self, _):
        self._update_kill_button_layout()

    def _on_process_change(self, _):
        self._model.update()
        self._update_state()
        self._update_kill_button_layout()
        self._update_clean_scratch_button_layout()

    def _on_kill_button_click(self, _):
        self._model.kill_process()
        self._update_kill_button_layout()

    def _on_update_results_button_click(self, _):
        self.node_view.node = None
        self.node_view.node = self._model.process_node

    def _on_clean_scratch_button_click(self, _):
        self._model.clean_remote_data()
        self._update_clean_scratch_button_layout()

    def _update_kill_button_layout(self):
        if not self.rendered:
            return
        if (
            not self._model.process
            or self._model.process_node.is_finished
            or self._model.process_node.is_excepted
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
        if self._model.process and self._model.process_node.is_terminated:
            self.clean_scratch_button.layout.display = "block"
        else:
            self.clean_scratch_button.layout.display = "none"

    def _update_state(self):
        if not self._model.process:
            self.state = self.State.INIT
        elif self._model.process_node.process_state in (
            ProcessState.CREATED,
            ProcessState.RUNNING,
            ProcessState.WAITING,
        ):
            self.state = self.State.ACTIVE
            self._model.process_info = PROCESS_RUNNING
        elif (
            self._model.process_node.process_state
            in (
                ProcessState.EXCEPTED,
                ProcessState.KILLED,
            )
            or self._model.process_node.is_failed
        ):
            self.state = self.State.FAIL
            self._model.process_info = PROCESS_EXCEPTED
        elif self._model.process_node.is_finished_ok:
            self.state = self.State.SUCCESS
            self._model.process_info = PROCESS_COMPLETED
