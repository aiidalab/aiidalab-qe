from __future__ import annotations

import ipywidgets as ipw

from aiida import orm
from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base.viewers import viewer as node_viewer

from .model import WorkChainStatusModel
from .tree import SimplifiedProcessTree, SimplifiedProcessTreeModel


class WorkChainStatusPanel(ResultsComponent[WorkChainStatusModel]):
    def __init__(self, model: WorkChainStatusModel, **kwargs):
        super().__init__(model=model, **kwargs)
        self.node_views = {}  # node-view cache
        self.node_view_loading_message = LoadingWidget("Loading node view")

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
        self.process_tree.observe(
            self._on_node_selection_change,
            "selected_nodes",
        )
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

        self.node_view_container = ipw.VBox()

        self.accordion = ipw.Accordion(
            children=[
                self.simplified_process_tree,
                ipw.VBox(
                    children=[
                        self.reset_button,
                        self.process_tree,
                        self.node_view_container,
                    ],
                ),
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

    def _on_monitor_counter_change(self, _):
        self._update_process_tree()

    def _on_accordion_change(self, change):
        if change["new"] == 0:
            self.simplified_process_tree.render()

    def _on_calculation_link_click(self, change):
        if selected_node_uuid := change["new"]:
            self.process_tree.value = selected_node_uuid
            self.accordion.selected_index = 1

    def _on_node_selection_change(self, change):
        self._update_node_view(change["new"])

    def _update_process_tree(self):
        if self.rendered:
            self.process_tree.update()

    def _update_node_view(self, nodes, refresh=False):
        """Update the node view based on the selected nodes.

        parameters
        ----------
        `nodes`: `list`
            List of selected nodes.
        `refresh`: `bool`, optional
            If True, the viewer will be refreshed.
            Occurs when user presses the "Update results" button.
        """

        if not nodes:
            return
        # only show the first selected node
        node = nodes[0]

        # check if the viewer is already added
        if node.uuid in self.node_views and not refresh:
            self.node_view = self.node_views[node.uuid]
        elif not isinstance(node, orm.WorkChainNode):
            self.node_view_container.children = [self.node_view_loading_message]
            self.node_view = node_viewer(node)
            self.node_views[node.uuid] = self.node_view
        else:
            self.node_view = ipw.HTML("No viewer available for this node.")

        self.node_view_container.children = [self.node_view]

    def _reset_process_tree(self, _):
        if not self.rendered:
            return
        self.process_tree.value = self._model.process_uuid
