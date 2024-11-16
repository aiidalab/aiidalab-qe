import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine import ProcessState
from aiidalab_widgets_base import (
    ProcessMonitor,
    ProcessNodesTreeWidget,
    WizardAppWidgetStep,
)
from aiidalab_widgets_base.viewers import viewer as node_viewer

from .model import ResultsStepModel
from .viewer import WorkChainViewer, WorkChainViewerModel

PROCESS_COMPLETED = "<h4>Workflow completed successfully!</h4>"
PROCESS_EXCEPTED = "<h4>Workflow is excepted!</h4>"
PROCESS_RUNNING = "<h4>Workflow is running!</h4>"


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppWidgetStep):
    def __init__(self, model: ResultsStepModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading results panel")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_process_change,
            "process_uuid",
        )

        self.rendered = False

        self.node_views = {}  # node-view cache

    def render(self):
        if self.rendered:
            return

        self.process_tree = ProcessNodesTreeWidget()
        self.process_tree.observe(
            self._on_node_selection_change,
            "selected_nodes",
        )
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.process_tree, "value"),
        )

        self.process_status = ipw.VBox(
            children=[
                self.process_tree,
            ],
        )

        self.process_monitor = ProcessMonitor(
            timeout=0.2,
            callbacks=[
                self.process_tree.update,
                self._update_status,
                self._update_state,
            ],
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

        # This starts the monitor on a separate thread
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
    def _on_state_change(self, _):
        self._update_kill_button_layout()

    def _on_process_change(self, _):
        self._model.update()
        self._update_state()
        self._update_kill_button_layout()
        self._update_clean_scratch_button_layout()

    def _on_node_selection_change(self, change):
        self._update_node_view(change["new"])

    def _on_kill_button_click(self, _):
        self._model.kill_process()
        self._update_kill_button_layout()

    def _on_update_results_button_click(self, _):
        self._update_node_view(self.process_tree.selected_nodes, refresh=True)

    def _on_clean_scratch_button_click(self, _):
        self._model.clean_remote_data()
        self._update_clean_scratch_button_layout()

    def _update_node_view(self, nodes, refresh=False):
        """Update the node view based on the selected nodes.

        parameters
        ----------
        `nodes`: `list`
            List of selected nodes.
        `refresh`: `bool`, optional
            If True, the viewer will be refreshed.
            Occurs when user presses the "Update results" button.
        """

        if not nodes:
            return
        # only show the first selected node
        node = nodes[0]

        # check if the viewer is already added
        if node.uuid in self.node_views and not refresh:
            self.node_view = self.node_views[node.uuid]
        elif not isinstance(node, orm.WorkChainNode):
            self.node_view = node_viewer(node)
            self.node_views[node.uuid] = self.node_view
        elif node.process_label == "QeAppWorkChain":
            self._create_workchain_viewer(node)

        self.process_status.children = [
            self.process_tree,
            self.node_view,
        ]

    def _create_workchain_viewer(self, node):
        model = WorkChainViewerModel()
        ipw.dlink(
            (self._model, "monitor_counter"),
            (model, "monitor_counter"),
        )
        self.node_view: WorkChainViewer = node_viewer(node, model=model)
        self.node_view.render()
        self.node_views[node.uuid] = self.node_view

    def _update_kill_button_layout(self):
        if not self.rendered:
            return
        process_node = self._model.fetch_process_node()
        if (
            not self._model.has_process
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
        process_node = self._model.fetch_process_node()
        if not self._model.has_process:
            self.state = self.State.INIT
        elif process_node.process_state in (
            ProcessState.CREATED,
            ProcessState.RUNNING,
            ProcessState.WAITING,
        ):
            self.state = self.State.ACTIVE
            self._model.process_info = PROCESS_RUNNING
        elif (
            process_node.process_state
            in (
                ProcessState.EXCEPTED,
                ProcessState.KILLED,
            )
            or process_node.is_failed
        ):
            self.state = self.State.FAIL
            self._model.process_info = PROCESS_EXCEPTED
        elif process_node.is_finished_ok:
            self.state = self.State.SUCCESS
            self._model.process_info = PROCESS_COMPLETED
        if self.state in (self.State.SUCCESS, self.State.FAIL):
            self._update_kill_button_layout()
            self._update_clean_scratch_button_layout()
