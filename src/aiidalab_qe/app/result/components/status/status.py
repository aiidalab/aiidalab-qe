import ipywidgets as ipw

from aiida import orm
from aiida.engine import ProcessState
from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.common.mixins import HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base.viewers import viewer as node_viewer

from .model import WorkChainStatusModel


class WorkChainStatusPanel(ResultsComponent[WorkChainStatusModel]):
    def __init__(self, model: WorkChainStatusModel, **kwargs):
        super().__init__(model=model, **kwargs)
        self.node_views = {}  # node-view cache
        self.node_view_loading_message = LoadingWidget("Loading node view")

    def _on_monitor_counter_change(self, _):
        self._update_simplified_view()
        self._update_process_tree()

    def _on_node_selection_change(self, change):
        self._update_node_view(change["new"])

    def _render(self):
        model = SimplifiedProcessTreeModel()
        self.simplified_process_tree = SimplifiedProcessTree(model=model)
        ipw.dlink(
            (self._model, "process_uuid"),
            (model, "process_uuid"),
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

        self.node_view_container = ipw.VBox()

        self.accordion = ipw.Accordion(
            children=[
                self.simplified_process_tree,
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

        self.accordion.observe(
            self._on_accordion_change,
            "selected_index",
        )

        self.children = [self.accordion]

    def _on_accordion_change(self, change):
        if change["new"] == 0:
            self.simplified_process_tree.render()

    def _update_simplified_view(self):
        if self.rendered:
            self.simplified_process_tree.update()

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


class TreeNode(ipw.VBox):
    def __init__(self, node, level=0, **kwargs):
        self.uuid = node.uuid
        self.level = level
        self.label = ipw.HTML(self._humanize_title(node))
        self.state = ""
        self.emoji = ipw.HTML()
        self.status = ipw.HTML()
        self.inspect = ipw.Button(
            description="Inspect",
            button_style="info",
            layout=ipw.Layout(width="fit-content", margin="0 0 0 5px"),
        )
        self.pks = set()
        self.title = ipw.HBox(
            children=[
                ipw.HTML(self._get_indentation(level)),
                self.emoji,
                self.label,
                self.status,
                self.inspect if isinstance(node, orm.CalcJobNode) else ipw.HTML(),
            ],
            layout=ipw.Layout(align_items="center"),
        )
        self.branches = ipw.VBox()
        super().__init__(
            children=[
                self.title,
                self.branches,
            ],
            **kwargs,
        )

    def update(self):
        node = orm.load_node(self.uuid)
        self._add_children(node)
        self.state = self._get_state(node)
        self.emoji.value = self._get_emoji(self.state)
        self.status.value = self._get_status(node)
        for branch in self.branches.children:
            if isinstance(branch, TreeNode):
                branch.update()

    def _add_children(self, node):
        for child in node.called:
            if child.pk in self.pks:
                continue
            if child.process_label == "BandsWorkChain":
                self._add_children(child)
            else:
                branch = TreeNode(child, level=self.level + 1)
                self.branches.children += (branch,)
                self.pks.add(child.pk)

    def _get_indentation(self, level=0):
        return "&nbsp;" * 8 * level

    def _get_emoji(self, state):
        return {
            "created": "🚀",
            "waiting": "💤",
            "running": "⏳",
            "finished": "✅",
            "killed": "💀",
            "excepted": "❌",
        }.get(state, "❓")

    def _get_status(self, node):
        return f"({self._get_tally(node)}{self.state})"

    def _get_tally(self, node):
        if not isinstance(node, orm.WorkflowNode):
            return ""
        inputs = node.get_metadata_inputs()
        processes = [key for key in inputs.keys() if key != "metadata"]
        total = len(processes)
        if node.process_label == "PwBaseWorkChain" and "kpoints" not in node.inputs:
            total += 1  # k-point grid generation
        if node.process_label == "PwBandsWorkChain":
            total += 1  # high-symmetry k-point generation
        finished = len(
            [
                child.process_state
                for child in node.called
                if child.process_state is ProcessState.FINISHED
            ]
        )
        return f"{finished}/{total} job{'s' if total > 1 else ''}; "

    def _get_state(self, node):
        if not hasattr(node, "process_state"):
            return "queued"
        state = node.process_state
        return (
            "running"
            if state is ProcessState.WAITING
            else state.value
            if state
            else "created"
        )

    def _humanize_title(self, node):
        if not hasattr(node, "process_label"):
            return "Unknown"
        title = node.process_label
        mappings = {
            "QeAppWorkChain": "Quantum ESPRESSO workflow",
            "BandsWorkChain": "Electronic band structure workflow",
            "PwBandsWorkChain": "Electronic band structure workflow",
            "PwRelaxWorkChain": "Structural relaxation workflow",
            "PwBaseWorkChain": "SCF workflow",
            "PhononWorkChain": "Phonons workflow",
            "PdosWorkChain": "Projected density of states workflow",
            "PwCalculation": "Run SCF cycle",
            "DosCalculation": "Compute density of states",
            "ProjwfcCalculation": "Compute projections",
            "create_kpoints_from_distance": "Generate K-points",
            "seekpath_structure_analysis": "Compute high-symmetry K-points",
        }
        return mappings.get(title, title)


class SimplifiedProcessTreeModel(Model, HasProcess):
    """"""


class SimplifiedProcessTree(ipw.VBox):
    def __init__(self, model: SimplifiedProcessTreeModel, **kwargs):
        super().__init__(**kwargs)
        self.add_class("simplified-process-tree")
        self._model = model
        self._model.observe(
            self._on_process_change,
            names="process_uuid",
        )
        self._model.observe(
            self._on_monitor_counter_change,
            "monitor_counter",
        )
        self.rendered = False

    def render(self):
        if self.rendered:
            return
        root = self._model.fetch_process_node()
        self.trunk = TreeNode(root)
        self.rendered = True
        self._update()
        self.children = [self.trunk]

    def _on_process_change(self, _):
        self._update()

    def _on_monitor_counter_change(self, _):
        self._update()

    def _update(self):
        if self.rendered:
            self.trunk.update()
