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

        self.kill_button = ipw.Button(
            description="Kill workchain",
            tooltip="Kill the below workchain.",
            button_style="danger",
            icon="window-close",
            layout=ipw.Layout(width="120px", height="40px"),
        )
        self.kill_button.on_click(self._on_click_kill_button)

        super().__init__([self.kill_button, self.process_status], **kwargs)

        self._update_kill_button_layout()

    def can_reset(self):
        "Do not allow reset while process is running."
        return self.state is not self.State.ACTIVE

    def reset(self):
        self.process = None

    def _update_state(self):
        """Based on the process state, update the state of the step."""
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
                self.kill_button.layout.display = "none"
            elif process.is_finished_ok:
                self.state = self.State.SUCCESS
                self.kill_button.layout.display = "none"

    def _update_kill_button_layout(self):
        """Update the layout of the kill button."""
        # If no process is selected, hide the button.
        if self.process is None or self.process == "":
            self.kill_button.layout.display = "none"
        else:
            process = orm.load_node(self.process)
            # If the process is finished or excepted, hide the button.
            if process.is_finished or process.is_excepted:
                self.kill_button.layout.display = "none"
            else:
                self.kill_button.layout.display = "block"

        # If the step is not activated, no point to click the button, so disable it.
        # Only enable it if the process is on (RUNNING, CREATED, WAITING).
        if self.state is self.State.ACTIVE:
            self.kill_button.disabled = False
        else:
            self.kill_button.disabled = True

    def _on_click_kill_button(self, _=None):
        """callback for the kill button.
        First kill the process, then update the kill button layout.
        """
        workchain = [orm.load_node(self.process)]
        control.kill_processes(workchain)

        # update the kill button layout
        self._update_kill_button_layout()

    @tl.observe("process")
    def _observe_process(self, _):
        """Callback for when the process is changed."""
        # The order of the following calls matters,
        # as the self.state is updated in the _update_state method.
        self._update_state()
        self._update_kill_button_layout()
