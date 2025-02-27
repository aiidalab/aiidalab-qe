from __future__ import annotations

import threading
import typing as t
from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.common.links import LinkType
from aiida.engine import ProcessState
from aiida.tools.graph.graph_traversers import traverse_graph
from aiidalab_qe.common.mixins import HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import LoadingWidget

from .state import STATE_ICONS


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
        if self.rendered:
            self._update()

    def _on_inspect(self, uuid: str):
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
        self.trunk.add_class("tree-trunk")
        self.trunk.initialize()
        self.trunk.expand()
        self._update()

        self.tree_container = ipw.VBox(
            children=[self.trunk],
        )
        self.tree_container.add_class("tree-container")

        self.children = [
            self.collapse_button,
            self.tree_container,
            ipw.HBox(
                children=[
                    ipw.HTML(
                        value="""
                        <div style="font-style: italic;">
                            *Actual number of jobs may exceed estimated total due to
                            error handling, dynamic sub-workflows, and/or other runtime
                            adjustments
                        </div>
                    """,
                    ),
                ],
                layout=ipw.Layout(
                    align_items="flex-start",
                    margin="5px 10px 0 0",
                ),
            ),
        ]

        self.rendered = True

    def _update(self):
        self.trunk.update()

    def _collapse_all(self, _):
        self.trunk.collapse(recursive=True)


ProcessNodeType = t.TypeVar("ProcessNodeType", bound=orm.ProcessNode)


class ProcessTreeNode(ipw.VBox, t.Generic[ProcessNodeType]):
    _TITLE_MAPPING = {
        "QeAppWorkChain": "Quantum ESPRESSO app workflow",
        "BandsWorkChain": "Electronic band structure workflow",
        "PwBandsWorkChain": "Electronic band structure workflow",
        "ProjwfcBandsWorkChain": "Electronic band structure workflow",
        "PwRelaxWorkChain": "Structure relaxation workflow manager",
        "PdosWorkChain": "Projected density of states workflow",
        "PwBaseWorkChain": {
            "scf": "SCF workflow",
            "nscf": "NSCF workflow",
            "bands": "Bands workflow",
            "relax": "Structure relaxation workflow",
            "md": "Molecular dynamics workflow",
        },
        "PwCalculation": {
            "scf": "Run SCF cycle",
            "nscf": "Run NSCF cycle",
            "bands": "Compute bands",
            "relax": "Optimize structure geometry",
            "md": "Run molecular dynamics simulation",
        },
        "DosCalculation": "Compute density of states",
        "ProjwfcCalculation": "Compute projections",
    }

    def __init__(
        self,
        node: ProcessNodeType,
        level: int = 0,
        on_inspect: t.Callable[[str], None] | None = None,
        **kwargs,
    ):
        self.uuid = node.uuid
        self.level = level
        self.on_inspect = on_inspect
        super().__init__(**kwargs)
        self._node: dict[int, ProcessNodeType] = {}  # thread_id: node

    @property
    def node(self) -> ProcessNodeType:
        if (tid := threading.get_ident()) not in self._node:
            node = orm.load_node(self.uuid)
            self._node[tid] = node
        return self._node[tid]

    def initialize(self):
        self._build_header()
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
        return STATE_ICONS.get(state, "❓")

    def _get_state(self):
        if not hasattr(self.node, "process_state"):
            return "queued"
        if self.node.is_failed:
            return "failed"
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
        if label in ("PwBaseWorkChain", "PwCalculation"):
            inputs = node.inputs.pw if label == "PwBaseWorkChain" else node.inputs
            calculation: str = inputs.parameters.get_dict()["CONTROL"]["calculation"]
            calculation = calculation.replace("vc-", "")
            return self._TITLE_MAPPING.get(label, {}).get(calculation, label)
        return self._TITLE_MAPPING.get(label, label)


class ProcessTreeBranches(ipw.VBox):
    def clear(self):
        self.children = []

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
    _adding_branches = False

    @property
    def metadata_inputs(self):
        # BACKWARDS COMPATIBILITY: originally this was added because "relax" was in the
        # metadata inputs regardless if running relaxation, as it was used also for the
        # summary to extract pw parameters. #1163 removes this dependency, thus allowing
        # the popping of the "relax" port if not running relaxation. However, this is
        # kept for backwards compatibility, for jobs ran before #1163.
        inputs = deepcopy(self.node.get_metadata_inputs()) or {}
        if "properties" in self.node.inputs:
            inputs = {
                key: value
                for key, value in inputs.items()
                if key in self.node.inputs.properties
            }
        return inputs

    @property
    def sub_processes(self):
        return set(self.metadata_inputs.keys()) - {"metadata"}

    @property
    def has_sub_workflows(self):
        return len(self.sub_processes) > 1

    @property
    def collapsed(self):
        return self.toggle.icon == "plus"

    def initialize(self):
        super().initialize()
        self.pks = set()
        self.branches = ProcessTreeBranches()
        self.branches.add_class("tree-node-branches")
        self.children += (self.branches,)
        self.expected_jobs = self._get_expected(self.metadata_inputs)
        self.update()

    def update(self):
        super().update()
        self.tally.value = self._get_tally()
        if not self.collapsed:
            self._add_branches()
        for branch in self.branches:
            branch.update()

    def clear(self):
        self.branches.clear()
        self.pks.clear()

    def expand(self, recursive=False):
        if self.collapsed:
            self.toggle.click()
        if recursive:
            for branch in self.branches:
                if isinstance(branch, WorkChainTreeNode):
                    branch.expand(recursive=True)

    def collapse(self, recursive: bool = False):
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

    def _add_branches(self):
        if self._adding_branches:
            return
        self._adding_branches = True
        self._add_branches_recursive()
        self._adding_branches = False

    def _add_branches_recursive(self, node: orm.ProcessNode | None = None):
        node = node or self.node
        for child in sorted(node.called, key=lambda child: child.ctime):
            if isinstance(child, orm.CalcFunctionNode):
                continue
            if child.pk in self.pks:
                continue
            if child.process_label in (
                "BandsWorkChain",
                "ProjwfcBaseWorkChain",
            ):
                self._add_branches_recursive(child)
            else:
                self._add_branch(child)

    def _add_branch(self, child: orm.ProcessNode):
        TreeNodeClass = (
            WorkChainTreeNode
            if isinstance(child, orm.WorkChainNode)
            else CalcJobTreeNode
        )
        branch = TreeNodeClass(
            node=child,
            level=self.level + 1,
            on_inspect=self.on_inspect,
        )
        branch.initialize()
        self.branches += branch
        self.pks.add(child.pk)

    def _get_tally(self):
        total = self.expected_jobs["count"]
        dynamic = self.expected_jobs["dynamic"]
        finished = self._count_finished()
        tally = f"{finished}/{total}"
        tally += "*" if dynamic else ""
        tally += " job" if total == 1 else " jobs"
        return tally

    def _count_finished(self):
        traverser = traverse_graph(
            starting_pks=[self.node.pk],
            links_forward=[
                LinkType.CALL_WORK,
                LinkType.CALL_CALC,
            ],
        )
        return (
            orm.QueryBuilder()
            .append(
                orm.CalcJobNode,
                filters=(
                    (orm.CalcJobNode.fields.pk.in_(traverser.get("nodes")))
                    & (orm.CalcJobNode.fields.process_state == "finished")
                    & (orm.CalcJobNode.fields.exit_status == 0)
                ),
            )
            .count()
        )

    def _get_expected(self, inputs: dict[str, dict]) -> dict:
        expected = {}
        count = 0
        dynamic = False

        for key, sub_inputs in inputs.items():
            if "metadata" in sub_inputs and "options" in sub_inputs["metadata"]:
                # This is a calculation
                count += 1
                dynamic = dynamic or key == "pw"
            elif key != "metadata":
                # This is a workflow
                nested = self._get_expected(sub_inputs)
                count += nested["count"]
                dynamic = dynamic or nested["dynamic"]
                expected[key] = nested

        expected |= {
            "count": count,
            "dynamic": dynamic,
        }

        return expected

    def _toggle_branches(self, _):
        if self.collapsed:
            self.branches.add_class("open")
            self.toggle.icon = "minus"
            self._add_branches()
        else:
            self.branches.remove_class("open")
            self.toggle.icon = "plus"


class CalcJobTreeNode(ProcessTreeNode[orm.CalcJobNode]):
    def initialize(self):
        super().initialize()
        self.update()

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
