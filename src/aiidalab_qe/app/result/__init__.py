import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine import ProcessState
from aiida.engine.processes import control
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
    process = tl.Unicode(allow_none=True)

    def __init__(self, model: ResultsModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading results panel")],
            **kwargs,
        )

        self._model = model

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
            icon="window-close",
            layout=ipw.Layout(width="120px", height="40px", display="none"),
        )
        ipw.dlink(
            (self, "state"),
            (self.kill_button, "disabled"),
            lambda state: state is not self.State.ACTIVE,
        )
        self.kill_button.on_click(self._on_click_kill_button)

        self.process_info = ipw.HTML()
        ipw.dlink(
            (self._model, "process_info"),
            (self.process_info, "value"),
        )

        self.children = [
            ipw.HBox(
                children=[
                    self.kill_button,
                    self.process_info,
                ]
            ),
            self.process_status,
        ]

        self._model.observe(
            self._on_process_change,
            "process",
        )

        ipw.dlink(
            (self._model, "process"),
            (self, "process"),
        )

        self.rendered = True

    def can_reset(self):
        return self.state is not self.State.ACTIVE

    def reset(self):
        self._model.reset()

    @tl.observe("state")
    def _on_state_change(self, _):
        self._update_kill_button_layout()

    def _on_process_change(self, _):
        self._update_state()
        self._update_kill_button_layout()

    def _on_click_kill_button(self, _=None):
        workchain = [orm.load_node(self.process)]
        control.kill_processes(workchain)
        self._update_kill_button_layout()

    def _update_kill_button_layout(self):
        if (
            self.process is None
            or self.process == ""
            or self.state
            in (
                self.State.SUCCESS,
                self.State.FAIL,
            )
        ):
            self.kill_button.layout.display = "none"
        else:
            process = orm.load_node(self.process)
            if process.is_finished or process.is_excepted:
                self.kill_button.layout.display = "none"
            else:
                self.kill_button.layout.display = "block"

    def _update_state(self):
        if self.process is None:
            self.state = self.State.INIT
        else:
            process = orm.load_node(self.process)
            process_state = process.process_state
            if process_state in (
                ProcessState.CREATED,
                ProcessState.RUNNING,
                ProcessState.WAITING,
            ):
                self.state = self.State.ACTIVE
                self._model.process_info = PROCESS_RUNNING
            elif (
                process_state
                in (
                    ProcessState.EXCEPTED,
                    ProcessState.KILLED,
                )
                or process.is_failed
            ):
                self.state = self.State.FAIL
                self._model.process_info = PROCESS_EXCEPTED
            elif process.is_finished_ok:
                self.state = self.State.SUCCESS
                self._model.process_info = PROCESS_COMPLETED
