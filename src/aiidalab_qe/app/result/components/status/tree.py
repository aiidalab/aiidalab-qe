from __future__ import annotations

import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine import BaseRestartWorkChain, ProcessState
from aiidalab_qe.common.mixins import HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import LoadingWidget


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
        if self._model.has_process:
            self._render()
        else:
            self.children = [ipw.HTML("Process node not yet available")]

    def _on_process_change(self, _):
        if not self.rendered:
            self._render()
        self._update()

    def _on_monitor_counter_change(self, _):
        self._update()

    def _on_inspect(self, uuid):
        self._model.clicked = None  # ensure event is triggered when label is reclicked
        self._model.clicked = uuid

    def _render(self):
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

    def _update(self):
        if self.rendered:
            self.trunk.update()

    def _collapse_all(self, _):
        self.trunk.collapse(recursive=True)


class ProcessTreeNode(ipw.VBox):
    _MAPPING = {
        "QeAppWorkChain": "Quantum ESPRESSO app workflow",
        "BandsWorkChain": "Electronic band structure workflow",
        "PwBandsWorkChain": "Electronic band structure workflow",
        "PwRelaxWorkChain": "Structural relaxation workflow",
        "PhononWorkChain": "Phonons workflow",
        "PdosWorkChain": "Projected density of states workflow",
        "DosCalculation": "Compute density of states",
        "ProjwfcCalculation": "Compute projections",
        "create_kpoints_from_distance": "Generate k-points",
        "seekpath_structure_analysis": "Compute high-symmetry k-points",
    }

    _PW_MAPPING = {
        "scf": {
            "PwBaseWorkChain": "SCF workflow",
            "PwCalculation": "Run SCF cycle",
        },
        "nscf": {
            "PwBaseWorkChain": "NSCF workflow",
            "PwCalculation": "Run NSCF cycle",
        },
        "bands": {
            "PwBaseWorkChain": "Bands workflow",
            "PwCalculation": "Compute bands",
        },
        "relax": {
            "PwBaseWorkChain": "Relaxation workflow",
            "PwCalculation": "Optimize structure geometry",
        },
    }

    def __init__(
        self,
        node: orm.ProcessNode,
        level: int = 0,
        on_inspect: t.Callable[[str], None] | None = None,
        **kwargs,
    ):
        if not (node and isinstance(node, orm.ProcessNode)):
            raise ValueError("Process node required")
        self.uuid = node.uuid
        self.on_inspect = on_inspect
        self._build_header(node, level)
        super().__init__(
            children=[self.header],
            **kwargs,
        )

    @property
    def process_node(self) -> orm.ProcessNode:
        return orm.load_node(self.uuid)  # type: ignore

    def update(self, node: orm.ProcessNode = None):
        node = node or self.process_node
        self.state.value = self._get_state(node)
        self.emoji.value = self._get_emoji(self.state.value)

    def _build_header(self, node: orm.ProcessNode, level: int):
        self.level = level
        self.title = self._get_human_readable_title(node)
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

    def _get_human_readable_title(self, node):
        if not hasattr(node, "process_label"):
            return "Unknown"
        label = node.process_label
        if label in ("PwBaseWorkChain", "PwCalculation"):
            inputs = node.inputs.pw if label == "PwBaseWorkChain" else node.inputs
            calculation: str = inputs.parameters.get_dict()["CONTROL"]["calculation"]
            mapping = self._PW_MAPPING[calculation]
        else:
            mapping = self._MAPPING
        return mapping.get(label, label)


class ProcessTreeBranches(ipw.VBox):
    def __len__(self):
        return len(self.children)

    def __getitem__(self, index: int) -> ProcessTreeNode:
        return self.children[index]  # type: ignore

    def __iter__(self) -> t.Iterator[ProcessTreeNode]:
        return iter(self.children)  # type: ignore

    def __iadd__(self, other: ProcessTreeNode):
        if not isinstance(other, ProcessTreeNode):
            raise TypeError("Right operand must be a TreeNode")
        self.children += (other,)
        return self


class WorkChainTreeNode(ProcessTreeNode):
    def __init__(
        self,
        node,
        level=0,
        on_inspect: t.Callable[[str], None] | None = None,
        **kwargs,
    ):
        super().__init__(node, level, on_inspect, **kwargs)
        self.pks = set()
        self.branches = ProcessTreeBranches()
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
        for branch in self.branches:
            branch.update()

    def expand(self, recursive=False):
        if self.collapsed:
            self.toggle.click()
        if recursive:
            for branch in self.branches:
                if isinstance(branch, WorkChainTreeNode):
                    branch.expand(recursive=True)

    def collapse(self, recursive=False):
        if not self.collapsed:
            self.toggle.click()
        if recursive:
            for branch in self.branches:
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
                self.branches += branch
                self.pks.add(child.pk)

    def _get_tally(self, node):
        total = self._get_total(node)
        finished = len(
            [
                child.process_state
                for child in node.called
                if child.process_state is ProcessState.FINISHED
            ]
        )
        tally = f"{finished}/{total}"
        tally += "*" if isinstance(node, BaseRestartWorkChain) else ""
        tally += " job" if total == 1 else " jobs"
        return tally

    def _get_total(self, node):
        if isinstance(node, BaseRestartWorkChain):
            total = len(node.called)
        else:
            inputs = node.get_metadata_inputs()
            processes = [key for key in inputs.keys() if key != "metadata"]
            total = len(processes)
        if node.process_label == "PwBaseWorkChain" and "kpoints" not in node.inputs:
            total += 1  # k-point grid generation
        if node.process_label == "PwBandsWorkChain":
            total += 1  # high-symmetry k-point generation
        return total

    def _toggle_branches(self, _):
        if self.collapsed:
            self.branches.add_class("open")
            self.toggle.icon = "minus"
            self._add_branches()
        else:
            self.branches.remove_class("open")
            self.toggle.icon = "plus"


class CalculationTreeNode(ProcessTreeNode):
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
