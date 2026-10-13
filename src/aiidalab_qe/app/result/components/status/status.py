from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base.viewers import AiidaNodeViewWidget

from .model import WorkChainStatusModel
from .paused import PausedProcessesModel, PausedProcessesTable
from .tree import SimplifiedProcessTree, SimplifiedProcessTreeModel


class WorkChainStatusPanel(ResultsComponent[WorkChainStatusModel]):
    def __init__(self, model: WorkChainStatusModel, **kwargs):
        super().__init__(model, **kwargs)

        model = SimplifiedProcessTreeModel()
        self.simplified_process_tree = SimplifiedProcessTree(model=model)
        self._model.add_model("tree", model)
        model.observe(
            self._on_calculation_link_click,
            "clicked",
        )

        self.paused_processes_model = PausedProcessesModel()
        self.paused_processes_table = PausedProcessesTable(
            model=self.paused_processes_model,
            on_inspect=self._on_paused_inspect,
        )
        self._model.add_model("paused", self.paused_processes_model)

    def _render(self):
        self.process_tree = ProcessNodesTreeWidget()
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.process_tree, "value"),
        )

        self.reset_button = ipw.Button(
            description="Reset to root",
            button_style="warning",
            icon="refresh",
            tooltip="Reseed the process tree with the root node",
            layout=ipw.Layout(width="fit-content"),
        )
        self.reset_button.on_click(self._reset_process_tree)

        self.node_view = AiidaNodeViewWidget()
        ipw.dlink(
            (self.process_tree, "selected_nodes"),
            (self.node_view, "node"),
            transform=lambda nodes: nodes[0] if nodes else None,
        )

        self.to_advanced_view_button = ipw.Button(
            description="View in advanced panel",
            button_style="primary",
            icon="eye",
            tooltip="Switch to the advanced view",
            layout=ipw.Layout(width="fit-content"),
        )
        self.to_advanced_view_button.on_click(self._switch_to_advanced_view)

        simplified_tree_section = ipw.Box(
            children=[
                ipw.VBox(
                    children=[
                        self.simplified_process_tree,
                    ],
                ),
                ipw.VBox(
                    children=[
                        self.to_advanced_view_button,
                        self.node_view,
                    ],
                ),
            ]
        )
        simplified_tree_section.add_class("simplified-view")

        advanced_tree_section = ipw.VBox(
            children=[
                self.reset_button,
                self.process_tree,
                self.node_view,
            ],
        )
        advanced_tree_section.add_class("advanced-view")

        paused_processes_section = ipw.VBox(
            children=[self.paused_processes_table],
        )

        self._sections: list[tuple[str, ipw.Widget, str]] = [
            (
                "paused",
                paused_processes_section,
                f"Paused Processes ({self.paused_processes_model.paused_count})",
            ),
            ("overview", simplified_tree_section, "Overview"),
            ("advanced", advanced_tree_section, "Advanced view"),
        ]
        self._section_indices = {
            name: index for index, (name, _, _) in enumerate(self._sections)
        }
        self.accordion = ipw.Accordion(
            children=[view for _, view, _ in self._sections],
            selected_index=None,
        )
        for name, _, title in self._sections:
            index = self._section_indices[name]
            self.accordion.set_title(index, title)

        self.accordion.observe(
            self._on_accordion_change,
            "selected_index",
        )
        self.paused_processes_model.observe(
            self._on_paused_count_change,
            "paused_count",
        )

        self.accordion.selected_index = self._section_indices["overview"]

        self.children = [self.accordion]

    def _post_render(self):
        self._select_tree_root()

    def _on_monitor_counter_change(self, _):
        if self.rendered:
            self.process_tree.update()

    def _on_accordion_change(self, change):
        if change["new"] == self._section_indices["overview"]:
            self.simplified_process_tree.render()

    def _on_calculation_link_click(self, change):
        if selected_node_uuid := change["new"]:
            self.process_tree.value = selected_node_uuid

    def _switch_to_advanced_view(self, _):
        self.accordion.selected_index = self._section_indices["advanced"]

    def _on_paused_count_change(self, change: dict):
        index = self._section_indices["paused"]
        self.accordion.set_title(index, f"Paused Processes ({change['new']})")

    def _on_paused_inspect(self, uuid):
        self.accordion.selected_index = self._section_indices["advanced"]
        self.process_tree.value = None
        self.process_tree.value = uuid

    def _select_tree_root(self):
        if self.rendered:
            self.process_tree.value = None
            self.process_tree.value = self._model.process_uuid

    def _reset_process_tree(self, _):
        if not self.rendered:
            return
        self.process_tree.value = self._model.process_uuid
