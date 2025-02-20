from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.common.process.tree import (
    SimplifiedProcessTree,
    SimplifiedProcessTreeModel,
)
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base.viewers import AiidaNodeViewWidget

from .model import WorkChainStatusModel


class WorkChainStatusPanel(ResultsComponent[WorkChainStatusModel]):
    def _render(self):
        model = SimplifiedProcessTreeModel()
        self.simplified_process_tree = SimplifiedProcessTree(model=model)
        ipw.dlink(
            (self._model, "process_uuid"),
            (model, "process_uuid"),
        )
        ipw.dlink(
            (self._model, "monitor_counter"),
            (model, "monitor_counter"),
        )
        model.observe(
            self._on_calculation_link_click,
            "clicked",
        )

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

        simplified_tree_container = ipw.VBox(
            children=[
                self.simplified_process_tree,
            ],
        )

        simplified_tree_node_view_container = ipw.VBox(
            children=[
                self.to_advanced_view_button,
                self.node_view,
            ],
        )

        simplified_view = ipw.Box(
            children=[
                simplified_tree_container,
                simplified_tree_node_view_container,
            ]
        )
        simplified_view.add_class("simplified-view")

        advanced_view = ipw.VBox(
            children=[
                self.reset_button,
                self.process_tree,
                self.node_view,
            ],
        )
        advanced_view.add_class("advanced-view")

        self.accordion = ipw.Accordion(
            children=[
                simplified_view,
                advanced_view,
            ],
            selected_index=None,
        )
        titles = [
            "Status overview",
            "Advanced status view",
        ]
        for i, title in enumerate(titles):
            self.accordion.set_title(i, title)

        self.accordion.observe(
            self._on_accordion_change,
            "selected_index",
        )

        self.accordion.selected_index = 0

        self.children = [self.accordion]

    def _post_render(self):
        self._select_tree_root()

    def _on_monitor_counter_change(self, _):
        if self.rendered:
            self.process_tree.update()

    def _on_accordion_change(self, change):
        if change["new"] == 0:
            self.simplified_process_tree.render()

    def _on_calculation_link_click(self, change):
        if selected_node_uuid := change["new"]:
            self.process_tree.value = selected_node_uuid

    def _switch_to_advanced_view(self, _):
        self.accordion.selected_index = 1

    def _select_tree_root(self):
        if self.rendered:
            self.process_tree.value = None
            self.process_tree.value = self._model.process_uuid

    def _reset_process_tree(self, _):
        if not self.rendered:
            return
        self.process_tree.value = self._model.process_uuid
