from __future__ import annotations

import typing as t

import ipywidgets as ipw

from aiida import orm
from aiida.engine import ProcessState
from aiidalab_qe.common.widgets import LoadingWidget

from .model import SimplifiedProcessTreeModel


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
        self.trunk.expand()
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
        self.trunk.collapse(recursive=True)


class TreeNode(ipw.VBox):
    _MAPPING = {
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

    @property
    def process_node(self):
        return orm.load_node(self.uuid)

    def update(self, node=None):
        node = node or self.process_node
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
        return self._MAPPING.get(title, title)


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
        node = node or self.process_node
        super().update(node)
        self.tally.value = self._get_tally(node)
        if not self.collapsed:
            self._add_branches(node)
        branch: TreeNode
        for branch in self.branches.children:
            branch.update()

    def expand(self, recursive=False):
        if self.collapsed:
            self.toggle.click()
        if recursive:
            for branch in self.branches.children:
                if isinstance(branch, WorkChainTreeNode):
                    branch.expand(recursive=True)

    def collapse(self, recursive=False):
        if not self.collapsed:
            self.toggle.click()
        if recursive:
            for branch in self.branches.children:
                if isinstance(branch, WorkChainTreeNode):
                    branch.collapse(recursive=True)

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

    def _add_branches(self, node=None):
        node = node or self.process_node
        for child in node.called:
            if child.pk in self.pks:
                continue
            if child.process_label == "BandsWorkChain":
                self._add_branches(child)
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
                branch.update()
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
            self._add_branches()
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
