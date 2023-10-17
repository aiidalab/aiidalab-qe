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
        
        self.kill_work_chains_button = ipw.Button(
            description="Kill workchain",
            tooltip="Kill the below workchain.",
            button_style="danger",
            icon="window-close",
            disabled=True,
            layout=ipw.Layout(width="120px",height="40px",display="none"),
        )
        self.kill_work_chains_button.on_click(self._on_click_kill_work_chain)
        self.kill_work_chains_box = ipw.HBox(
            children=[self.kill_work_chains_button],
            layout=ipw.Layout(justify_content="flex-end")
            )

        super().__init__([self.kill_work_chains_box,self.process_status], **kwargs)

    def can_reset(self):
        "Do not allow reset while process is running."
        return self.state is not self.State.ACTIVE

    def reset(self):
        self.process = None

    def _update_state(self):
        if self.process is None:
            self.state = self.State.INIT
            self.kill_work_chains_button.display = "none"
        else:
            self.kill_work_chains_button.layout.display = "block"
            process = orm.load_node(self.process)
            process_state = process.process_state
            if process_state in (
                ProcessState.CREATED,
                ProcessState.RUNNING,
                ProcessState.WAITING,
            ):
                self.state = self.State.ACTIVE
                self.kill_work_chains_button.disabled = False
            elif (
                process_state in (ProcessState.EXCEPTED, ProcessState.KILLED)
                or process.is_failed
            ):
                self.state = self.State.FAIL
                self.kill_work_chains_button.disabled = True
            elif process.is_finished_ok:
                self.state = self.State.SUCCESS
                self.kill_work_chains_button.disabled = True
                
    def _on_click_kill_work_chain(self, _=None):
        workchain = [orm.load_node(self.process)]
        control.kill_processes(workchain)
        self._update_state()

    @tl.observe("process")
    def _observe_process(self, change):
        self._update_state()