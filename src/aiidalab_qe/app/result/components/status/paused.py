from __future__ import annotations

import typing as t
from dataclasses import dataclass

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine.processes import control
from aiidalab_qe.app.result.utils import ResultsSubModel


@dataclass(frozen=True, slots=True)
class PausedProcess:
    uuid: str
    pk: int
    label: str
    status: str


PausedProcessRowType = tuple[ipw.HBox, ipw.HTML, ipw.HTML, ipw.Button]

PK_LAYOUT = ipw.Layout(width="60px")
LABEL_LAYOUT = ipw.Layout(flex="1 1 0px", min_width="0")
REASON_LAYOUT = ipw.Layout(flex="2 1 0px", min_width="0")
ACTIONS_LAYOUT = ipw.Layout(width="72px")
ACTION_BUTTON_LAYOUT = ipw.Layout(width="36px")
FULL_WIDTH_LAYOUT = ipw.Layout(width="100%")
PLAY_ALL_LAYOUT = ipw.Layout(width="fit-content", margin="0 0 0 auto")
TABLE_LAYOUT = ipw.Layout(grid_gap="8px")


class PausedProcessesModel(ResultsSubModel):
    """Tracks paused processes (calculations and workflows) within the monitored process."""

    paused_count = tl.Int(0)
    error_message = tl.Unicode("")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.paused_processes: tuple[PausedProcess, ...] = ()
        self.observe(lambda _: self.update(), "monitor_counter")

    def play(self, uuid: str):
        """Attempt to resume a single paused process."""
        try:
            control.play_processes([orm.load_node(uuid)])
        except Exception as exception:
            self.error_message = str(exception)
        else:
            self.error_message = ""
        self.update()

    def play_all(self):
        """Attempt to resume all currently paused processes."""
        try:
            processes = [
                orm.load_node(process.uuid) for process in self.paused_processes
            ]
            control.play_processes(processes)
        except Exception as exception:
            self.error_message = str(exception)
        else:
            self.error_message = ""
        self.update()

    def reset(self):
        self.paused_processes = ()
        self.paused_count = 0
        self.error_message = ""
        self.daemon_is_running = False
        self.daemon_status_known = False

    def _update(self, specific=""):
        if not self.has_process:
            self.paused_processes = ()
        else:
            self.paused_processes = tuple(
                PausedProcess(
                    uuid=node.uuid,
                    pk=node.pk,
                    label=node.label or node.process_label,
                    status=node.process_status or "Paused",
                )
                for node in (self.process, *self.process.called_descendants)
                if (
                    isinstance(node, orm.ProcessNode)
                    and node.paused
                    and not node.is_terminated
                )
            )
        self.paused_count = len(self.paused_processes)


class PausedProcessesTable(ipw.VBox):
    """Table of paused processes, with actions to inspect or resume them."""

    def __init__(
        self,
        model: PausedProcessesModel,
        on_inspect: t.Callable[[str], None] | None = None,
        **kwargs,
    ):
        self._model = model
        self.on_inspect = on_inspect

        self.table = ipw.VBox()
        self.alert = ipw.HTML()
        self.play_all_button = ipw.Button(
            description="Play all",
            icon="play",
            button_style="success",
            tooltip="Resume all paused processes in this workflow",
            layout=PLAY_ALL_LAYOUT,
        )
        self.play_all_button.on_click(lambda _: self._model.play_all())

        self._header = self._build_header_row()
        self._empty_message = ipw.HTML("<b>No paused processes</b>")
        self._rows: dict[str, PausedProcessRowType] = {}
        self._row_states: dict[str, tuple[str, str]] = {}
        self._row_order: tuple[str, ...] | None = None

        super().__init__(
            children=[
                self.play_all_button,
                self.table,
                self.alert,
            ],
            layout=TABLE_LAYOUT,
            **kwargs,
        )

        self._model.observe(
            self._on_monitor_counter_change,
            "monitor_counter",
        )
        self._model.observe(
            self._on_paused_count_change,
            "paused_count",
        )
        self._model.observe(
            self._on_error_message_change,
            "error_message",
        )
        self._model.observe(
            self._on_daemon_status_change,
            [
                "daemon_is_running",
                "daemon_status_known",
            ],
        )

        self._render_table()
        self._update_play_all_button()

    def _on_monitor_counter_change(self, _):
        self._render_table()
        self._update_play_all_button()

    def _on_paused_count_change(self, _):
        self._render_table()
        self._update_play_all_button()

    def _on_error_message_change(self, change):
        message = change["new"]
        self.alert.value = (
            f"""
            <div class="alert alert-danger">
                {message}
            </div>
            """
            if message
            else ""
        )

    def _on_daemon_status_change(self, _):
        self._update_play_buttons()
        self._update_play_all_button()

    def _render_table(self):
        processes = {process.uuid: process for process in self._model.paused_processes}
        self._sync_rows(processes)
        self._set_table_rows(tuple(processes))

    def _sync_rows(self, processes: dict[str, PausedProcess]):
        """Synchronize the internal row representations with the given nodes.

        This method ensures that the internal dictionaries `_rows` and `_row_states`
        are consistent with the provided `processes` dictionary. Rows corresponding to
        removed nodes are deleted, and new rows are created for newly added nodes.
        Existing rows are updated if their state has changed.
        """
        for uuid in self._rows.keys() - processes.keys():
            self._rows.pop(uuid)
            self._row_states.pop(uuid)

        for uuid, process in processes.items():
            state = self._get_row_state(process)
            if row := self._rows.get(uuid):
                if self._row_states[uuid] != state:
                    self._update_row(row, state)
            else:
                self._rows[uuid] = self._build_row(process)
            self._row_states[uuid] = state

    def _get_row_state(self, process: PausedProcess) -> tuple[str, str]:
        return process.label, process.status

    def _update_row(
        self,
        row: PausedProcessRowType,
        state: tuple[str, str],
    ):
        _, label, reason, _ = row
        label.value = state[0]
        reason.value = state[1]

    def _set_table_rows(self, order: tuple[str, ...]):
        """Set the table rows in the specified order.

        If the order has not changed, this method does nothing.
        Otherwise, it updates the table's children to reflect the new order,
        including the header and any empty message if there are no rows.
        """
        if order == self._row_order:
            return
        self._row_order = order
        self.table.children = (
            (self._empty_message,)
            if not order
            else (self._header, *(self._rows[uuid][0] for uuid in order))
        )
        self._update_play_buttons()
        self._update_play_all_button()

    def _build_header_row(self):
        return ipw.HBox(
            children=[
                ipw.HTML("<b>PK</b>", layout=PK_LAYOUT),
                ipw.HTML("<b>Label</b>", layout=LABEL_LAYOUT),
                ipw.HTML("<b>Reason</b>", layout=REASON_LAYOUT),
                ipw.HTML("<b>Actions</b>", layout=ACTIONS_LAYOUT),
            ],
            layout=FULL_WIDTH_LAYOUT,
        )

    def _build_row(self, process: PausedProcess) -> PausedProcessRowType:
        label = ipw.HTML(process.label, layout=LABEL_LAYOUT)
        reason = ipw.HTML(process.status, layout=REASON_LAYOUT)

        goto_button = ipw.Button(
            icon="share",
            tooltip="Go to process in advanced status view",
            layout=ACTION_BUTTON_LAYOUT,
        )
        goto_button.on_click(lambda _, uuid=process.uuid: self._on_goto_click(uuid))

        play_button = ipw.Button(
            icon="play",
            button_style="success",
            tooltip=(
                "Resume the paused process"
                if self._model.daemon_status_known and self._model.daemon_is_running
                else "Resume is unavailable until the AiiDA daemon status is known"
            ),
            disabled=(
                not self._model.daemon_status_known or not self._model.daemon_is_running
            ),
            layout=ACTION_BUTTON_LAYOUT,
        )
        play_button.on_click(lambda _, uuid=process.uuid: self._model.play(uuid))

        row = ipw.HBox(
            children=[
                ipw.HTML(f"<b>{process.pk}</b>", layout=ipw.Layout(width="60px")),
                label,
                reason,
                goto_button,
                play_button,
            ],
            layout=FULL_WIDTH_LAYOUT,
        )
        return row, label, reason, play_button

    def _update_play_buttons(self):
        for _, _, _, play_button in self._rows.values():
            play_button.disabled = (
                not self._model.daemon_status_known or not self._model.daemon_is_running
            )
            play_button.tooltip = (
                "Resume the paused process"
                if self._model.daemon_status_known and self._model.daemon_is_running
                else "Resume is unavailable until the AiiDA daemon status is known"
            )

    def _update_play_all_button(self):
        self.play_all_button.layout.display = (
            "block" if self._model.paused_processes else "none"
        )
        self.play_all_button.disabled = (
            not self._model.paused_processes
            or not self._model.daemon_status_known
            or not self._model.daemon_is_running
        )
        self.play_all_button.tooltip = (
            "Resume all paused processes in this workflow"
            if self._model.daemon_status_known and self._model.daemon_is_running
            else "Resume is unavailable until the AiiDA daemon status is known"
        )

    def _on_goto_click(self, uuid: str):
        if self.on_inspect:
            self.on_inspect(uuid)
