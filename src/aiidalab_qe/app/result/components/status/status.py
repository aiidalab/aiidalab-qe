from __future__ import annotations

import typing as t

import ipywidgets as ipw
import traitlets as tl

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


class TreeNode(ipw.VBox):
    def __init__(
        self,
        node,
        level=0,
        on_inspect: t.Callable[[str], None] | None = None,
        **kwargs,
    ):
        self.uuid = node.uuid
        self.on_inspect = on_inspect
        self._build_header(node, level)
        super().__init__(
            children=[self.header],
            **kwargs,
        )

    def update(self, node=None):
        node = node or orm.load_node(self.uuid)
        self.state.value = self._get_state(node)
        self.emoji.value = self._get_emoji(self.state.value)

    def _build_header(self, node, level):
        self.level = level
        self.title = self._humanize_title(node)
        self.indentation = self._get_indentation(level)
        self.emoji = ipw.HTML()
        self.state = ipw.HTML()
        self.header = ipw.HBox()
        self.header.add_class("tree-node-header")

    def _get_indentation(self, level=0):
        return ipw.HTML(layout=ipw.Layout(width=f"{22 * level}px"))

    def _get_emoji(self, state):
        return {
            "created": "🚀",
            "waiting": "💤",
            "running": "⏳",
            "finished": "✅",
            "killed": "💀",
            "excepted": "❌",
        }.get(state, "❓")

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


class WorkChainTreeNode(TreeNode):
    def __init__(
        self,
        node,
        level=0,
        on_inspect: t.Callable[[str], None] | None = None,
        **kwargs,
    ):
        super().__init__(node, level, on_inspect, **kwargs)
        self.pks = set()
        self.branches = ipw.VBox()
        self.branches.add_class("tree-node-branches")
        self.children += (self.branches,)

    @property
    def collapsed(self):
        return self.toggle.icon == "plus"

    def update(self, node=None):
        node = node or orm.load_node(self.uuid)
        super().update(node)
        self.tally.value = self._get_tally(node)
        self._add_children(node)
        branch: TreeNode
        for branch in self.branches.children:
            branch.update()

    def collapse(self):
        if not self.collapsed:
            self.toggle.click()
        for branch in self.branches.children:
            if isinstance(branch, WorkChainTreeNode):
                branch.collapse()

    def _build_header(self, node, level):
        super()._build_header(node, level)
        self.toggle = ipw.Button(icon="plus")
        self.toggle.add_class("tree-node-toggle")
        self.toggle.on_click(self._toggle_branches)
        self.label = ipw.HTML(self.title)
        self.tally = ipw.HTML()
        self.header.children = [
            self.indentation,
            self.toggle,
            self.emoji,
            self.label,
            ipw.HTML(" | "),
            self.state,
            ipw.HTML(" | "),
            self.tally,
        ]

    def _add_children(self, node):
        for child in node.called:
            if child.pk in self.pks:
                continue
            if child.process_label == "BandsWorkChain":
                self._add_children(child)
            else:
                TreeNodeClass = (
                    WorkChainTreeNode
                    if isinstance(child, orm.WorkflowNode)
                    else CalculationTreeNode
                )
                branch = TreeNodeClass(
                    child,
                    level=self.level + 1,
                    on_inspect=self.on_inspect,
                )
                self.branches.children += (branch,)
                self.pks.add(child.pk)

    def _get_tally(self, node):
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

    def _toggle_branches(self, _):
        if self.collapsed:
            self.branches.add_class("open")
            self.toggle.icon = "minus"
        else:
            self.branches.remove_class("open")
            self.toggle.icon = "plus"


class CalculationTreeNode(TreeNode):
    def _build_header(self, node, level):
        super()._build_header(node, level)
        self.label = ipw.Button(
            description=self.title,
            tooltip="click to view calculation",
        )
        self.label.add_class("calculation-link")
        self.label.on_click(self._on_label_click)
        self.header.children = [
            self.indentation,
            self.emoji,
            self.label,
            ipw.HTML(" | "),
            self.state,
        ]

    def _on_label_click(self, _):
        if self.on_inspect is None:
            return
        self.on_inspect(self.uuid)


class SimplifiedProcessTreeModel(Model, HasProcess):
    clicked = tl.Unicode(None, allow_none=True)


class SimplifiedProcessTree(ipw.VBox):
    def __init__(self, model: SimplifiedProcessTreeModel, **kwargs):
        self.loading_message = LoadingWidget("Loading process tree")
        super().__init__(
            children=[self.loading_message],
            **kwargs,
        )
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
        self.collapse_button = ipw.Button(
            description="Collapse all",
            button_style="warning",
            icon="minus",
            tooltip="Collapse all branches",
            layout=ipw.Layout(
                width="fit-content",
                margin="0 0 10px 0",
            ),
        )
        self.collapse_button.on_click(self._collapse_all)
        root = self._model.fetch_process_node()
        self.trunk = WorkChainTreeNode(node=root, on_inspect=self._on_inspect)
        self.rendered = True
        self._update()
        self.children = [
            self.collapse_button,
            self.trunk,
        ]

    def _on_process_change(self, _):
        self._update()

    def _on_monitor_counter_change(self, _):
        self._update()

    def _on_inspect(self, uuid):
        self._model.clicked = None  # ensure event is triggered
        self._model.clicked = uuid

    def _update(self):
        if self.rendered:
            self.trunk.update()

    def _collapse_all(self, _):
        self.trunk.collapse()
