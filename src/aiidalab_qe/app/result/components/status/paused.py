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

    def _on_error_message_change(self, change):
        message = change["new"]
        self.alert.value = (
            f'<div class="alert alert-danger">{message}</div>' if message else ""
        )

    def _on_daemon_status_change(self, _):
        self._render_daemon_warning()
        self._render_table()

    def _render_daemon_warning(self):
        self.daemon_warning.value = (
            '<div class="alert alert-warning">'
            "⚠️ The AiiDA daemon is not running. Process resuming is disabled."
            "</div>"
            if self._model.paused_count and not self._model.daemon_is_running
            else ""
        )

    def _render_table(self):
        if not self._model.nodes:
            self.table.children = [ipw.HTML("<b>No paused processes</b>")]
            return
        self.table.children = [
            self._build_header_row(),
            *[self._build_row(node) for node in self._model.nodes],
        ]

    def _build_header_row(self):
        return ipw.HBox(
            children=[
                ipw.HTML("<b>PK</b>", layout=ipw.Layout(width="60px")),
                ipw.HTML("<b>Label</b>", layout=ipw.Layout(width="200px")),
                ipw.HTML("<b>Reason</b>", layout=ipw.Layout(flex="1")),
                ipw.HTML("<b>Actions</b>", layout=ipw.Layout(width="72px")),
            ],
        )

    def _build_row(self, node: orm.ProcessNode):
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

        return ipw.HBox(
            children=[
                ipw.HTML(f"<b>{node.pk}</b>", layout=ipw.Layout(width="60px")),
                ipw.HTML(
                    node.label or node.process_label,
                    layout=ipw.Layout(width="200px"),
                ),
                ipw.HTML(
                    node.process_status or "Paused",
                    layout=ipw.Layout(flex="1"),
                ),
                goto_button,
                play_button,
            ],
        )

    def _on_goto_click(self, uuid: str):
        if self.on_inspect:
            self.on_inspect(uuid)
