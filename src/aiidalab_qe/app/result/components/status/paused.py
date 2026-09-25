from __future__ import annotations

import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine.processes import control
from aiidalab_qe.app.result.utils import ResultsSubModel


class PausedProcessesModel(ResultsSubModel):
    """Tracks paused processes (calculations and workflows) within the monitored process."""

    paused_count = tl.Int(0)
    error_message = tl.Unicode("")

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

    def _update(self, specific=""):
        self.nodes = (
            [
                node
                for node in (self.process, *self.process.called_descendants)
                if isinstance(node, orm.ProcessNode) and node.paused
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
        self.alert = ipw.HTML()

        super().__init__(children=[self.table, self.alert], **kwargs)

        self._model.observe(
            self._on_monitor_counter_change,
            "monitor_counter",
        )
        self._model.observe(
            self._on_error_message_change,
            "error_message",
        )

        self._render_table()

    def _on_monitor_counter_change(self, _):
        self._render_table()

    def _on_error_message_change(self, change):
        message = change["new"]
        self.alert.value = (
            f'<div class="alert alert-danger">{message}</div>' if message else ""
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
            tooltip="Replay the paused process",
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
