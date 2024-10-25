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

PROCESS_COMPLETED = "<h4>Workflow completed successfully!</h4>"
PROCESS_EXCEPTED = "<h4>Workflow is excepted!</h4>"
PROCESS_RUNNING = "<h4>Workflow is running!</h4>"


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
            icon="stop",
            layout=ipw.Layout(width="120px", display="none", margin="0px 20px 0px 0px"),
        )
        self.kill_button.on_click(self._on_click_kill_button)

        self.clean_scratch_button = ipw.Button(
            description="Clean remote data",
            tooltip="Clean the remote folders of the workchain.",
            button_style="danger",
            icon="trash",
            layout=ipw.Layout(width="150px", display="none", margin="0px 20px 0px 0px"),
        )
        self.clean_scratch_button.on_click(self._on_click_clean_scratch_button)
        self.update_result_button = ipw.Button(
            description="Update results tabs",
            tooltip="Trigger the update of the results tabs.",
            button_style="success",
            icon="refresh",
            layout=ipw.Layout(
                width="150px", display="block", margin="0px 20px 0px 0px"
            ),
        )
        self.update_result_button.on_click(self._on_click_update_result_button)

        self.process_info = ipw.HTML()

        super().__init__(
            [
                self.process_info,
                ipw.HBox(
                    children=[
                        self.kill_button,
                        self.update_result_button,
                        self.clean_scratch_button,
                    ]
                ),
                self.process_status,
            ],
            **kwargs,
        )

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
                self.process_info.value = PROCESS_RUNNING
            elif (
                process_state in (ProcessState.EXCEPTED, ProcessState.KILLED)
                or process.is_failed
            ):
                self.state = self.State.FAIL
                self.process_info.value = PROCESS_EXCEPTED
            elif process.is_finished_ok:
                self.state = self.State.SUCCESS
                self.process_info.value = PROCESS_COMPLETED
            # trigger the update of kill and clean button.
            if self.state in [self.State.SUCCESS, self.State.FAIL]:
                self._update_kill_button_layout()
                self._update_clean_scratch_button_layout()

    def _update_kill_button_layout(self):
        """Update the layout of the kill button."""
        # If no process is selected, hide the button.
        if self.process is None or self.process == "":
            self.kill_button.layout.display = "none"
        else:
            process = orm.load_node(self.process)
            # If the process is terminated, hide the button.
            if process.is_terminated:
                self.kill_button.layout.display = "none"
            else:
                self.kill_button.layout.display = "block"

        # If the step is not activated, no point to click the button, so disable it.
        # Only enable it if the process is on (RUNNING, CREATED, WAITING).
        if self.state is self.State.ACTIVE:
            self.kill_button.disabled = False
        else:
            self.kill_button.disabled = True

    def _update_clean_scratch_button_layout(self):
        """Update the layout of the kill button."""
        # The button is hidden by default, but if we load a new process, we hide again.
        if not self.process:
            self.clean_scratch_button.layout.display = "none"
        else:
            process = orm.load_node(self.process)
            # If the process is terminated, show the button.
            if process.is_terminated:
                self.clean_scratch_button.layout.display = "block"
            else:
                self.clean_scratch_button.layout.display = "none"

            # If the scratch is already empty, we should deactivate the button.
            # not sure about the performance if descendants are several.
            cleaned_bool = []
            for called_descendant in process.called_descendants:
                if isinstance(called_descendant, orm.CalcJobNode):
                    try:
                        cleaned_bool.append(
                            called_descendant.outputs.remote_folder.is_empty
                        )
                    except Exception:
                        pass
            self.clean_scratch_button.disabled = all(cleaned_bool)

    def _on_click_kill_button(self, _=None):
        """callback for the kill button.
        First kill the process, then update the kill button layout.
        """
        workchain = [orm.load_node(self.process)]
        control.kill_processes(workchain)

        # update the kill button layout
        self._update_kill_button_layout()

    def _on_click_clean_scratch_button(self, _=None):
        """callback for the clean scratch button.
        First clean the remote folders, then update the clean button layout.
        """
        process = orm.load_node(self.process)

        for called_descendant in process.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()
                except Exception:
                    pass

        # update the kill button layout
        self._update_clean_scratch_button_layout()

    def _on_click_update_result_button(self, _=None):
        """Trigger the update of the results tabs."""
        # change the node to trigger the update of the view.
        self.node_view.node = None
        self.node_view.node = orm.load_node(self.process)

    @tl.observe("process")
    def _observe_process(self, _):
        """Callback for when the process is changed."""
        # The order of the following calls matters,
        # as the self.state is updated in the _update_state method.
        self._update_state()
        self._update_kill_button_layout()
        self._update_clean_scratch_button_layout()
