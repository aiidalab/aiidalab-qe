from __future__ import annotations

import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine.daemon.client import get_daemon_client
from aiida.engine.processes import control
from aiidalab_qe.app.result.utils import ResultsSubModel


class PausedProcessesModel(ResultsSubModel):
    """Tracks paused processes (calculations and workflows) within the monitored process."""

    paused_count = tl.Int(0)
    error_message = tl.Unicode("")
    daemon_is_running = tl.Bool(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nodes: list[orm.ProcessNode] = []
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

    def reset(self):
        self.process_uuid = None
        self.nodes = []
        self.paused_count = 0
        self.error_message = ""
        self.daemon_is_running = False

    def _update(self, specific=""):
        try:
            self.daemon_is_running = get_daemon_client().is_daemon_running
        except Exception:
            self.daemon_is_running = False
        self.nodes = (
            [
                node
                for node in (self.process, *self.process.called_descendants)
                if (
                    isinstance(node, orm.ProcessNode)
                    and node.paused
                    and not node.is_terminated
                )
            ]
            if self.has_process
            else []
        )
        self.paused_count = len(self.nodes)


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
        self.daemon_warning = ipw.HTML()
        self.alert = ipw.HTML()

        self._header = self._build_header_row()
        self._empty_message = ipw.HTML("<b>No paused processes</b>")
        self._rows: dict[str, tuple[ipw.HBox, ipw.HTML, ipw.HTML, ipw.Button]] = {}
        self._row_states: dict[str, tuple[str, str]] = {}
        self._row_order: tuple[str, ...] | None = None

        super().__init__(
            children=[
                self.daemon_warning,
                self.table,
                self.alert,
            ],
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
            "daemon_is_running",
        )

        self._render_table()
        self._render_daemon_warning()

    def _on_monitor_counter_change(self, _):
        self._render_table()
        self._render_daemon_warning()

    def _on_paused_count_change(self, _):
        self._render_table()
        self._render_daemon_warning()

    def _on_error_message_change(self, change):
        message = change["new"]
        self.alert.value = (
            f'<div class="alert alert-danger">{message}</div>' if message else ""
        )

    def _on_daemon_status_change(self, _):
        self._render_daemon_warning()
        self._update_play_buttons()

    def _render_daemon_warning(self):
        self.daemon_warning.value = (
            '<div class="alert alert-warning">'
            "⚠️ The AiiDA daemon is not running. Process resuming is disabled."
            "</div>"
            if self._model.paused_count and not self._model.daemon_is_running
            else ""
        )

    def _render_table(self):
        nodes = {node.uuid: node for node in self._model.nodes}
        self._sync_rows(nodes)
        self._set_table_rows(tuple(nodes))

    def _sync_rows(self, nodes: dict[str, orm.ProcessNode]):
        """Synchronize the internal row representations with the given nodes.

        This method ensures that the internal dictionaries `_rows` and `_row_states`
        are consistent with the provided `nodes` dictionary. Rows corresponding to
        removed nodes are deleted, and new rows are created for newly added nodes.
        Existing rows are updated if their state has changed.
        """
        for uuid in self._rows.keys() - nodes.keys():
            self._rows.pop(uuid)
            self._row_states.pop(uuid)

        for uuid, node in nodes.items():
            state = self._get_row_state(node)
            if row := self._rows.get(uuid):
                if self._row_states[uuid] != state:
                    self._update_row(row, state)
            else:
                self._rows[uuid] = self._build_row(node)
            self._row_states[uuid] = state

    def _get_row_state(self, node: orm.ProcessNode) -> tuple[str, str]:
        return (
            node.label or node.process_label,
            node.process_status or "Paused",
        )

    def _update_row(
        self,
        row: tuple[ipw.HBox, ipw.HTML, ipw.HTML, ipw.Button],
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

    def _build_header_row(self):
        return ipw.HBox(
            children=[
                ipw.HTML("<b>PK</b>", layout=ipw.Layout(width="60px")),
                ipw.HTML("<b>Label</b>", layout=ipw.Layout(width="200px")),
                ipw.HTML("<b>Reason</b>", layout=ipw.Layout(flex="1")),
                ipw.HTML("<b>Actions</b>", layout=ipw.Layout(width="72px")),
            ],
        )

    def _build_row(
        self,
        node: orm.ProcessNode,
    ) -> tuple[ipw.HBox, ipw.HTML, ipw.HTML, ipw.Button]:
        label = ipw.HTML(
            node.label or node.process_label,
            layout=ipw.Layout(width="200px"),
        )
        reason = ipw.HTML(
            node.process_status or "Paused",
            layout=ipw.Layout(flex="1"),
        )

        goto_button = ipw.Button(
            icon="share",
            tooltip="Go to process in advanced status view",
            layout=ipw.Layout(width="36px"),
        )
        goto_button.on_click(lambda _, uuid=node.uuid: self._on_goto_click(uuid))

        play_button = ipw.Button(
            icon="play",
            button_style="success",
            tooltip=(
                "Resume the paused process"
                if self._model.daemon_is_running
                else "Resume is unavailable because the AiiDA daemon is not running"
            ),
            disabled=not self._model.daemon_is_running,
            layout=ipw.Layout(width="36px"),
        )
        play_button.on_click(lambda _, uuid=node.uuid: self._model.play(uuid))

        row = ipw.HBox(
            children=[
                ipw.HTML(f"<b>{node.pk}</b>", layout=ipw.Layout(width="60px")),
                label,
                reason,
                goto_button,
                play_button,
            ],
        )
        return row, label, reason, play_button

    def _update_play_buttons(self):
        for _, _, _, play_button in self._rows.values():
            play_button.disabled = not self._model.daemon_is_running
            play_button.tooltip = (
                "Resume the paused process"
                if self._model.daemon_is_running
                else "Resume is unavailable because the AiiDA daemon is not running"
            )

    def _on_goto_click(self, uuid: str):
        if self.on_inspect:
            self.on_inspect(uuid)
