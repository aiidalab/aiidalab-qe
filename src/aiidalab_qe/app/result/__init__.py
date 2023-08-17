import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiida.engine import ProcessState
from aiidalab_widgets_base import (
    AiidaNodeViewWidget,
    ProcessMonitor,
    ProcessNodesTreeWidget,
    WizardAppWidgetStep,
)

# trigger registration of the viewer widget:
from .workchain_viewer import WorkChainViewer  # noqa: F401


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppWidgetStep):
    process = tl.Unicode(allow_none=True)

    def __init__(self, **kwargs):
        self.process_tree = ProcessNodesTreeWidget()
        ipw.dlink(
            (self, "process"),
            (self.process_tree, "value"),
        )

        self.node_view = AiidaNodeViewWidget(layout={"width": "auto", "height": "auto"})
        ipw.dlink(
            (self.process_tree, "selected_nodes"),
            (self.node_view, "node"),
            transform=lambda nodes: nodes[0] if nodes else None,
        )
        self.process_status = ipw.VBox(children=[self.process_tree, self.node_view])

        # Setup process monitor
        self.process_monitor = ProcessMonitor(
            timeout=0.2,
            callbacks=[
                self.process_tree.update,
                self._update_state,
            ],
        )
        ipw.dlink((self, "process"), (self.process_monitor, "value"))

        super().__init__([self.process_status], **kwargs)

    def can_reset(self):
        "Do not allow reset while process is running."
        return self.state is not self.State.ACTIVE

    def reset(self):
        self.process = None

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
            elif (
                process_state in (ProcessState.EXCEPTED, ProcessState.KILLED)
                or process.is_failed
            ):
                self.state = self.State.FAIL
            elif process.is_finished_ok:
                self.state = self.State.SUCCESS

    @tl.observe("process")
    def _observe_process(self, change):
        self._update_state()
