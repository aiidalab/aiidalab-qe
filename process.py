"""Widgets for the monitoring of processes."""
import itertools
import os
from dataclasses import dataclass
from threading import Thread, Event, Lock

import traitlets
import ipywidgets as ipw
from IPython.display import clear_output
from IPython.display import display
from aiida.cmdline.utils.ascii_vis import calc_info
from aiida.cmdline.utils.query.calculation import CalculationQueryBuilder
from aiida.engine import ProcessState
from aiida.orm import CalcFunctionNode
from aiida.orm import CalcJobNode
from aiida.orm import Node
from aiida.orm import ProcessNode
from aiida.orm import WorkChainNode
from aiida.orm import load_node
from ipytree import Node as TreeNode
from ipytree import Tree


def get_calc_job_output(process):
    from aiidalab_widgets_base.process import get_running_calcs

    previous_calc_id = None
    num_lines = 0

    while not process.is_sealed:
        calc = None
        for calc in get_running_calcs(process):
            if calc.id == previous_calc_id:
                break
        else:
            if calc:
                previous_calc_id = calc.id

        if calc and "remote_folder" in calc.outputs:
            f_path = os.path.join(
                calc.outputs.remote_folder.get_remote_path(),
                calc.attributes["output_filename"],
            )
            if os.path.exists(f_path):
                with open(f_path) as fobj:
                    new_lines = fobj.readlines()[num_lines:]
                    num_lines += len(new_lines)
                    yield from new_lines


class ProgressBarWidget(ipw.VBox):
    """A bar showing the proggress of a process."""

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        self.correspondance = {
            None: (0, "warning"),
            "created": (0, "info"),
            "running": (1, "info"),
            "waiting": (1, "info"),
            "killed": (2, "danger"),
            "excepted": (2, "danger"),
            "finished": (2, "success"),
        }
        self.bar = ipw.IntProgress(  # pylint: disable=blacklisted-name
            value=0,
            min=0,
            max=2,
            step=1,
            bar_style="warning",  # 'success', 'info', 'warning', 'danger' or ''
            orientation="horizontal",
            layout=ipw.Layout(width="auto"),
        )
        self.state = ipw.HTML(
            description="Calculation state:",
            value="",
            style={"description_width": "100px"},
        )
        super().__init__(children=[self.state, self.bar], **kwargs)

    @traitlets.observe("process")
    def update(self, _=None):
        """Update the bar."""
        self.bar.value, self.bar.bar_style = self.correspondance[self.current_state]
        if self.current_state is None:
            self.state.value = "N/A"
        else:
            self.state.value = self.current_state.capitalize()

    @property
    def current_state(self):
        if self.process is not None:
            return self.process.process_state.value


class ProcessMonitor(traitlets.HasTraits):
    """Monitor a process and execute callback functions at specified intervals."""

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, callbacks=None, timeout=None, **kwargs):
        self.callbacks = [] if callbacks is None else list(callbacks)
        self.timeout = 0.1 if timeout is None else timeout

        self._monitor_thread = None
        self._monitor_thread_stop = Event()
        self._monitor_thread_lock = Lock()

        super().__init__(**kwargs)

    @traitlets.observe("process")
    def _observe_process(self, change):
        process = change["new"]
        if process is None or process.id != getattr(change["old"], "id", None):
            with self.hold_trait_notifications():
                with self._monitor_thread_lock:
                    # stop thread
                    if self._monitor_thread is not None:
                        self._monitor_thread_stop.set()
                        self._monitor_thread.join()

                    # reset output
                    for callback, _ in self.callbacks:
                        callback(None)

                    # start monitor thread
                    self._monitor_thread_stop.clear()
                    process_id = getattr(process, "id", None)
                    self._monitor_thread = Thread(
                        target=self._monitor_process, args=(process_id,)
                    )
                    self._monitor_thread.start()

    def _monitor_process(self, process_id):
        self._monitor_thread_stop.wait(
            timeout=10 * self.timeout
        )  # brief delay to increase app stability

        process = None if process_id is None else load_node(process_id)

        iterations = itertools.count()
        while not (process is None or process.is_sealed):

            iteration = next(iterations)
            for callback, period in self.callbacks:
                if iteration % period == 0:
                    callback(process_id)

            if self._monitor_thread_stop.wait(timeout=self.timeout):
                break

        # Final update:
        for callback, _ in self.callbacks:
            callback(process_id)


class WorkChainSelector(ipw.HBox):

    # The PK of a 'aiida.workflows:quantumespresso.pw.bands' WorkChainNode.
    value = traitlets.Int(allow_none=True)

    # When this trait is set to a positive value, the work chains are automatically
    # refreshed every `auto_refresh_interval` seconds.
    auto_refresh_interval = traitlets.Int()  # seconds

    # Indicate whether the widget is currently updating the work chain options.
    busy = traitlets.Bool(read_only=True)

    # Note: We use this class as a singleton to reset the work chains selector
    # widget to its default stage (no work chain selected), because we cannot
    # use `None` as setting the widget's value to None will lead to "no selection".
    _NO_PROCESS = object()

    FMT_WORKCHAIN = "{wc.pk:6}{wc.ctime:>10}\t{wc.state:<16}\t{wc.formula}"

    def __init__(self, **kwargs):
        self.work_chains_selector = ipw.Dropdown(
            description="WorkChain",
            options=[("New calculation...", self._NO_PROCESS)],
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
        )
        ipw.dlink(
            (self.work_chains_selector, "value"),
            (self, "value"),
            transform=lambda pk: None if pk is self._NO_PROCESS else pk,
        )

        self.refresh_work_chains_button = ipw.Button(description="Refresh")
        self.refresh_work_chains_button.on_click(self.refresh_work_chains)

        self._refresh_lock = Lock()
        self._refresh_thread = None
        self._stop_refresh_thread = Event()
        self._update_auto_refresh_thread_state()

        super().__init__(
            children=[self.work_chains_selector, self.refresh_work_chains_button],
            **kwargs,
        )

    @dataclass
    class WorkChainData:
        pk: int
        ctime: str
        state: str
        formula: str

    @classmethod
    def find_work_chains(cls):
        builder = CalculationQueryBuilder()
        filters = builder.get_filters(
            process_label="PwBandsWorkChain",
        )
        query_set = builder.get_query_set(
            filters=filters,
            order_by={"ctime": "desc"},
        )
        projected = builder.get_projected(
            query_set, projections=["pk", "ctime", "state"]
        )

        for process in projected[1:]:
            pk = process[0]
            formula = load_node(pk).inputs.structure.get_formula()
            yield cls.WorkChainData(formula=formula, *process)

    @traitlets.default("busy")
    def _default_busy(self):
        return True

    @traitlets.observe("busy")
    def _observe_busy(self, change):
        for child in self.children:
            child.disabled = change["new"]

    def refresh_work_chains(self, _=None):
        with self._refresh_lock:
            try:
                self.set_trait("busy", True)  # disables the widget

                with self.hold_trait_notifications():
                    # We need to restore the original value, because it may be reset due to this issue:
                    # https://github.com/jupyter-widgets/ipywidgets/issues/2230
                    original_value = self.work_chains_selector.value

                    self.work_chains_selector.options = [
                        ("New calculation...", self._NO_PROCESS)
                    ] + [
                        (self.FMT_WORKCHAIN.format(wc=wc), wc.pk)
                        for wc in self.find_work_chains()
                    ]

                    self.work_chains_selector.value = original_value
            finally:
                self.set_trait("busy", False)  # reenable the widget

    def _auto_refresh_loop(self):
        while True:
            self.refresh_work_chains()
            if self._stop_refresh_thread.wait(timeout=self.auto_refresh_interval):
                break

    def _update_auto_refresh_thread_state(self):
        if self.auto_refresh_interval > 0 and self._refresh_thread is None:
            # start thread
            self._stop_refresh_thread.clear()
            self._refresh_thread = Thread(target=self._auto_refresh_loop)
            self._refresh_thread.start()

        elif self.auto_refresh_interval <= 0 and self._refresh_thread is not None:
            # stop thread
            self._stop_refresh_thread.set()
            self._refresh_thread.join(timeout=30)
            self._refresh_thread = None

    @traitlets.default("auto_refresh_interval")
    def _default_auto_refresh_interval(self):
        return 10  # seconds

    @traitlets.observe("auto_refresh_interval")
    def _observe_auto_refresh_interval(self, change):
        if change["new"] != change["old"]:
            self._update_auto_refresh_thread_state()

    @traitlets.observe("value")
    def _observe_value(self, change):
        if change["old"] == change["new"]:
            return

        new = self._NO_PROCESS if change["new"] is None else change["new"]

        if new not in {pk for _, pk in self.work_chains_selector.options}:
            self.refresh_work_chains()

        self.work_chains_selector.value = new


class AiidaNodeTreeNode(TreeNode):
    def __init__(self, pk, name, **kwargs):
        self.pk = pk
        self.nodes_registry = dict()
        super().__init__(name=name, **kwargs)

    @traitlets.default("opened")
    def _default_openend(self):
        return False


class AiidaProcessNodeTreeNode(AiidaNodeTreeNode):
    def __init__(self, pk, **kwargs):
        self.outputs_node = AiidaOutputsTreeNode(name="outputs", parent_pk=pk)
        super().__init__(pk=pk, **kwargs)


class WorkChainProcessTreeNode(AiidaProcessNodeTreeNode):
    icon = traitlets.Unicode("chain").tag(sync=True)


class CalcJobTreeNode(AiidaProcessNodeTreeNode):
    icon = traitlets.Unicode("gears").tag(sync=True)


class CalcFunctionTreeNode(AiidaProcessNodeTreeNode):
    icon = traitlets.Unicode("gear").tag(sync=True)


class AiidaOutputsTreeNode(TreeNode):
    icon = traitlets.Unicode("folder").tag(sync=True)
    disabled = traitlets.Bool(True).tag(sync=True)

    def __init__(self, name, parent_pk, **kwargs):
        self.parent_pk = parent_pk
        self.nodes_registry = dict()
        super().__init__(name=name, **kwargs)


class UnknownTypeTreeNode(AiidaNodeTreeNode):
    icon = traitlets.Unicode("file").tag(sync=True)


class NodesTreeWidget(ipw.Output):
    """A tree widget for the structured representation of a nodes graph."""

    nodes = traitlets.Tuple().tag(trait=traitlets.Instance(Node))
    selected_nodes = traitlets.Tuple(read_only=True).tag(trait=traitlets.Instance(Node))

    PROCESS_STATE_STYLE = {
        ProcessState.EXCEPTED: "danger",
        ProcessState.FINISHED: "success",
        ProcessState.KILLED: "warning",
        ProcessState.RUNNING: "info",
        ProcessState.WAITING: "info",
    }

    PROCESS_STATE_STYLE_DEFAULT = "default"

    NODE_TYPE = {
        WorkChainNode: WorkChainProcessTreeNode,
        CalcFunctionNode: CalcFunctionTreeNode,
        CalcJobNode: CalcJobTreeNode,
    }

    def __init__(self, **kwargs):
        self._tree = Tree()
        self._tree.observe(self._observe_tree_selected_nodes, ["selected_nodes"])

        super().__init__(**kwargs)

    def _refresh_output(self):
        # There appears to be a bug in the ipytree implementation that sometimes
        # causes the output to not be properly cleared. We therefore refresh the
        # displayed tree upon change of the process trait.
        with self:
            clear_output()
            display(self._tree)

    def _observe_tree_selected_nodes(self, change):
        return self.set_trait(
            "selected_nodes",
            tuple(
                load_node(pk=node.pk) for node in change["new"] if hasattr(node, "pk")
            ),
        )

    def _convert_to_tree_nodes(self, old_nodes, new_nodes):
        "Convert nodes into tree nodes while re-using already converted nodes."
        old_nodes_ = {node.pk: node for node in old_nodes}
        assert len(old_nodes_) == len(old_nodes)  # no duplicated nodes

        for node in new_nodes:
            if node.pk in old_nodes_:
                yield old_nodes_[node.pk]
            else:
                yield self._to_tree_node(node, opened=True)

    @traitlets.observe("nodes")
    def _observe_nodes(self, change):
        self._tree.nodes = list(
            sorted(
                self._convert_to_tree_nodes(
                    old_nodes=self._tree.nodes, new_nodes=change["new"]
                ),
                key=lambda node: node.pk,
            )
        )
        self.update()
        self._refresh_output()

    @classmethod
    def _to_tree_node(cls, node, name=None, **kwargs):
        """Convert an AiiDA node to a tree node."""
        if name is None:
            if isinstance(node, ProcessNode):
                name = calc_info(node)
            else:
                name = str(node)
        return cls.NODE_TYPE.get(type(node), UnknownTypeTreeNode)(
            pk=node.pk, name=name, **kwargs
        )

    @classmethod
    def _find_called(cls, root):
        assert isinstance(root, AiidaProcessNodeTreeNode)
        process_node = load_node(root.pk)
        called = process_node.called
        called.sort(key=lambda p: p.ctime)
        for node in called:
            if node.pk not in root.nodes_registry:
                root.nodes_registry[node.pk] = cls._to_tree_node(
                    node, name=calc_info(node)
                )
            yield root.nodes_registry[node.pk]

    @classmethod
    def _find_outputs(cls, root):
        assert isinstance(root, AiidaOutputsTreeNode)
        process_node = load_node(root.parent_pk)
        outputs = {k: process_node.outputs[k] for k in process_node.outputs}
        for key in sorted(outputs.keys(), key=lambda k: outputs[k].pk):
            output_node = outputs[key]
            if output_node.pk not in root.nodes_registry:
                root.nodes_registry[output_node.pk] = cls._to_tree_node(
                    output_node, name=f"{key}<{output_node.pk}>"
                )
            yield root.nodes_registry[output_node.pk]

    @classmethod
    def _find_children(cls, root):
        """Find all children of the provided AiiDA node."""
        if isinstance(root, AiidaProcessNodeTreeNode):
            yield root.outputs_node
            yield from cls._find_called(root)
        elif isinstance(root, AiidaOutputsTreeNode):
            yield from cls._find_outputs(root)

    @classmethod
    def _build_tree(cls, root):
        """Recursively build a tree nodes graph for a given tree node."""
        root.nodes = [cls._build_tree(child) for child in cls._find_children(root)]
        return root

    @classmethod
    def _walk_tree(cls, root):
        """Breadth-first search of the node tree."""
        yield root
        for node in root.nodes:
            yield from cls._walk_tree(node)

    def _update_tree_node(self, tree_node):
        if isinstance(tree_node, AiidaProcessNodeTreeNode):
            process_node = load_node(tree_node.pk)
            tree_node.name = calc_info(process_node)
            tree_node.icon_style = self.PROCESS_STATE_STYLE.get(
                process_node.process_state, self.PROCESS_STATE_STYLE_DEFAULT
            )

    def update(self, _=None):
        """Refresh nodes based on the latest state of the root process and its children."""
        for root_node in self._tree.nodes:
            self._build_tree(root_node)
            for tree_node in self._walk_tree(root_node):
                self._update_tree_node(tree_node)

    def find_node(self, pk):
        for node in self._walk_tree(self._tree):
            if getattr(node, "pk", None) == pk:
                return node
        raise KeyError(pk)


class ProcessNodesTreeWidget(ipw.VBox):
    """A tree widget for the structured representation of a process graph.

    Args:
        refresh_period:
            The time period in between updates to the process tree view in seconds.
    """

    process = traitlets.Instance(ProcessNode, allow_none=True)
    selected_nodes = traitlets.Tuple(read_only=True).tag(trait=traitlets.Instance(Node))

    def __init__(self, refresh_period=0.2, **kwargs):
        self._tree = NodesTreeWidget()
        self._tree.observe(self._observe_tree_selected_nodes, ["selected_nodes"])

        if refresh_period > 0:
            self._process_monitor = ProcessMonitor(
                timeout=refresh_period,
                callbacks=[
                    (self.update, 1),
                ],
            )
            ipw.dlink((self, "process"), (self._process_monitor, "process"))
        else:
            self._process_monitor = None  # externally managed

        super().__init__(children=[self._tree], **kwargs)

    def _observe_tree_selected_nodes(self, change):
        self.set_trait("selected_nodes", change["new"])

    def update(self, _=None):
        self._tree.update()

    @traitlets.observe("process")
    def _observe_process(self, change):
        process = change["new"]
        if process:
            self._tree.nodes = [process]
            self._tree.find_node(process.pk).selected = True
        else:
            self._tree.nodes = []
