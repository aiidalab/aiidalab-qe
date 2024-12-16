import ipywidgets as ipw

from aiida import orm
from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base.viewers import viewer as node_viewer

from .model import WorkChainStatusModel


class WorkChainStatusPanel(ResultsComponent[WorkChainStatusModel]):
    def __init__(self, model: WorkChainStatusModel, **kwargs):
        super().__init__(model=model, **kwargs)
        self.node_views = {}  # node-view cache
        self.node_view_loading_message = LoadingWidget("Loading node view")

    def _render(self):
        self.simplified_status_view = ipw.HTML("Coming soon")

        self.process_tree = ProcessNodesTreeWidget()
        self.process_tree.observe(
            self._on_node_selection_change,
            "selected_nodes",
        )
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.process_tree, "value"),
        )

        self.node_view_container = ipw.VBox()

        self.accordion = ipw.Accordion(
            children=[
                self.simplified_status_view,
                ipw.VBox(
                    children=[
                        self.process_tree,
                        self.node_view_container,
                    ],
                ),
            ],
            selected_index=1,
        )
        titles = [
            "Status overview",
            "Advanced status view",
        ]
        for i, title in enumerate(titles):
            self.accordion.set_title(i, title)

        self.children = [self.accordion]

    def _on_monitor_counter_change(self, _):
        self._update_process_tree()

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
