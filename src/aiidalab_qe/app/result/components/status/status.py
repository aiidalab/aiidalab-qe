import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine import ProcessState
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

    def _on_monitor_counter_change(self, _):
        self._update_simplified_view()
        self._update_process_tree()

    def _on_node_selection_change(self, change):
        self._update_node_view(change["new"])

    def _render(self):
        self.simplified_process_tree = SimplifiedProcessTreeWidget()
        ipw.dlink(
            (self._model, "process_uuid"),
            (self.simplified_process_tree, "process_uuid"),
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

        self.children = [self.accordion]

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


class SimplifiedProcessTreeWidget(ipw.HTML):
    process_uuid = tl.Unicode(None, allow_none=True)

    def __init__(self, value=None, **kwargs):
        super().__init__(value, **kwargs)
        self.observe(
            self._on_process_change,
            "process_uuid",
        )

    def update(self):
        root = orm.load_node(self.process_uuid)
        simplified_tree = self._build_node(root)
        self.value = self._format(simplified_tree)

    def _on_process_change(self, _):
        self.update()

    def _build_node(self, node):
        tree_node = {
            "children": [],
        }

        for child in list(node.called):
            child_node = self._build_node(child)
            if child.process_label == "BandsWorkChain" and child_node["children"]:
                tree_node["children"].append(child_node["children"][0])
            else:
                tree_node["children"].append(child_node)

        tree_node["title"] = self._prepare_title(node)

        return tree_node

    def _prepare_title(self, node):
        progress = (
            self._prepare_progress_string(node)
            if isinstance(node, orm.WorkflowNode)
            else ""
        )
        title = self._humanize_title(node)
        state = self._get_state(node)
        status = f"({f'{progress}; ' if progress else ''}{state})"
        emoji = {
            "created": "🚀",
            "waiting": "💤",
            "running": "⏳",
            "finished": "✅",
            "killed": "💀",
            "excepted": "❌",
        }.get(state, "❓")
        return f"{emoji} {title} {status}"

    def _prepare_progress_string(self, node):
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
        return f"{finished}/{total} job{'s' if total > 1 else ''}"

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

    def _format(self, tree, level=0):
        indent = "&nbsp;" * 8 * level
        output = f"{indent}{tree['title']}"
        for child in tree["children"]:
            output += f"<br>{self._format(child, level + 1)}"
        return output
