from __future__ import annotations

import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine import ProcessState
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
        self.trunk.initialize()
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


ProcessNodeType = t.TypeVar("ProcessNodeType", bound=orm.ProcessNode)


class ProcessTreeNode(ipw.VBox, t.Generic[ProcessNodeType]):
    _MAPPING = {
        "QeAppWorkChain": "Quantum ESPRESSO app workflow",
        "PwBandsWorkChain": "Electronic band structure workflow",
        "PwRelaxWorkChain": "Structural relaxation workflow",
        "PdosWorkChain": "Projected density of states workflow",
        "DosCalculation": "Compute density of states",
        "ProjwfcCalculation": "Compute projections",
    }

    _PW_MAPPING = {
        "scf": "Run SCF cycle",
        "nscf": "Run NSCF cycle",
        "bands": "Compute bands",
        "relax": "Optimize structure geometry",
    }

    def __init__(
        self,
        node: ProcessNodeType,
        level: int = 0,
        on_inspect: t.Callable[[str], None] | None = None,
        **kwargs,
    ):
        # if not (node and isinstance(node, orm.ProcessNode)):
        #     raise ValueError("Process node required")
        self.uuid = node.uuid
        self.level = level
        self.on_inspect = on_inspect
        super().__init__(**kwargs)

    @property
    def node(self) -> ProcessNodeType:
        return orm.load_node(self.uuid)  # type: ignore

    def initialize(self):
        self._build_header()
        self.update()
        self.children = [self.header]

    def update(self):
        self.state.value = self._get_state()
        self.emoji.value = self._get_emoji(self.state.value)

    def _build_header(self):
        self.title = self._get_human_readable_title()
        self.indentation = self._get_indentation()
        self.emoji = ipw.HTML()
        self.state = ipw.HTML()
        self.header = ipw.HBox()
        self.header.add_class("tree-node-header")

    def _get_indentation(self):
        return ipw.HTML(layout=ipw.Layout(width=f"{22 * self.level}px"))

    def _get_emoji(self, state):
        return {
            "created": "🚀",
            "waiting": "💤",
            "running": "⏳",
            "finished": "✅",
            "killed": "💀",
            "excepted": "❌",
        }.get(state, "❓")

    def _get_state(self):
        if not hasattr(self.node, "process_state"):
            return "queued"
        state = self.node.process_state
        return (
            "running"
            if state is ProcessState.WAITING
            else state.value
            if state
            else "created"
        )

    def _get_human_readable_title(self):
        node = self.node
        if not (label := node.process_label):
            return "Unknown"
        if label not in ("PwBaseWorkChain", "PwCalculation"):
            return self._MAPPING.get(label, label)
        inputs = node.inputs.pw if label == "PwBaseWorkChain" else node.inputs
        calculation: str = inputs.parameters.get_dict()["CONTROL"]["calculation"]
        return self._PW_MAPPING[calculation]


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


class WorkChainTreeNode(ProcessTreeNode[orm.WorkChainNode]):
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

    @property
    def metadata_inputs(self):
        return self.node.get_metadata_inputs()

    @property
    def has_sub_workflows(self):
        return len(set(self.metadata_inputs.keys()) - {"metadata"}) > 1

    @property
    def collapsed(self):
        return self.toggle.icon == "plus"

    def initialize(self):
        super().initialize()
        self.children += (self.branches,)

    def update(self):
        super().update()
        self.tally.value = self._get_tally()
        if not self.collapsed:
            self._add_branches()
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

    def _build_header(self):
        super()._build_header()
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
        node = node or self.node
        for child in node.called:
            if child.pk in self.pks or isinstance(child, orm.CalcFunctionNode):
                continue
            TreeNodeClass = (
                WorkChainTreeNode
                if isinstance(child, orm.WorkChainNode)
                else CalcJobTreeNode
            )
            branch = TreeNodeClass(
                child,
                level=self.level + 1,
                on_inspect=self.on_inspect,
            )
            if isinstance(branch, WorkChainTreeNode) and not branch.has_sub_workflows:
                self._add_branches(child)
            else:
                branch.initialize()
                self.branches += branch
                self.pks.add(child.pk)

    def _get_tally(self):
        inputs = self.metadata_inputs
        total, dynamic = self._count_calculations(inputs)
        finished = self._count_finished_calculations(self.node)
        tally = f"{finished}/{total}"
        tally += "*" if dynamic and total != finished else ""
        tally += " job" if total == 1 else " jobs"
        return tally

    def _count_finished_calculations(self, node):
        count = 0
        for child in node.called:
            if isinstance(child, orm.WorkChainNode):
                count += self._count_finished_calculations(child)
            elif (
                isinstance(child, orm.CalcJobNode)
                and child.process_state is ProcessState.FINISHED
            ):
                count += 1
        return count

    def _count_calculations(
        self,
        inputs: dict[str, dict],
        parent: str | None = None,
    ) -> tuple[int, bool]:
        """Count all calculations in the metadata inputs dictionary of a workflow.

        Recursively count occurrences of the "metadata" key in the metadata inputs
        dictionary. Also check if any "metadata" key has a "pw" parent and mark the
        workflow as dynamic if found. The function ends branch exploration once
        "metadata" is encountered.

        Parameters
        ----------
        `inputs`: `dict[str, dict]`
            The metadata inputs dictionary of a workflow.
        `parent`: `str`, optional
            The parent key of the current metadata inputs dictionary.

        Returns
        -------
        `tuple[int, book]`
            The number of calculations found and if the workflow is dynamic.
        """
        count = 0
        dynamic = False

        for key, sub_inputs in inputs.items():
            if key == "metadata" and "options" in sub_inputs:  # a calculation
                count += 1
                dynamic = parent == "pw"
                return count, dynamic
            if isinstance(sub_inputs, dict):
                sub_count, sub_dynamic = self._count_calculations(sub_inputs, key)
                count += sub_count
                dynamic = dynamic or sub_dynamic

        return count, dynamic

    def _toggle_branches(self, _):
        if self.collapsed:
            self.branches.add_class("open")
            self.toggle.icon = "minus"
            self._add_branches()
        else:
            self.branches.remove_class("open")
            self.toggle.icon = "plus"


class CalcJobTreeNode(ProcessTreeNode[orm.CalcJobNode]):
    def _build_header(self):
        super()._build_header()
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
