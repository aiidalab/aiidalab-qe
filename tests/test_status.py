from __future__ import annotations

import threading
import time
from unittest.mock import Mock

import ipywidgets as ipw
import pytest

from aiida import orm
from aiida.common.links import LinkType
from aiida.engine import ProcessState
from aiida.engine.processes import control
from aiidalab_qe.app.result.components.status import (
    WorkChainStatusModel,
    WorkChainStatusPanel,
)
from aiidalab_qe.app.result.components.status.paused import (
    PausedProcessesModel,
    PausedProcessesTable,
)
from aiidalab_qe.app.result.components.status.tree import (
    TITLE_MAPPING,
    CalcJobTreeNode,
    SimplifiedProcessTree,
    SimplifiedProcessTreeModel,
    WorkChainTreeNode,
)


def mock_calcjob(label):
    calcjob = orm.CalcJobNode()
    calcjob.set_process_state(ProcessState.FINISHED)
    calcjob.set_exit_status(0)
    calcjob.set_process_label(label)
    return calcjob


def mock_workchain(label):
    workchain = orm.WorkChainNode()
    workchain.set_process_state(ProcessState.FINISHED)
    workchain.set_exit_status(0)
    workchain.set_process_label(label)
    return workchain


def create_process_graph():
    root = orm.WorkChainNode()
    root.set_process_label("QeAppWorkChain")
    root.set_process_state(ProcessState.RUNNING)

    child = orm.WorkChainNode()
    child.set_process_label("PwBaseWorkChain")
    child.set_process_state(ProcessState.RUNNING)

    calculation = orm.CalcJobNode()
    calculation.set_process_label("PwCalculation")
    calculation.set_process_state(ProcessState.RUNNING)

    child.base.links.add_incoming(
        root,
        link_type=LinkType.CALL_WORK,
        link_label="child",
    )
    calculation.base.links.add_incoming(
        child,
        link_type=LinkType.CALL_CALC,
        link_label="calculation",
    )

    root.store()
    child.store()
    calculation.store()
    return root, child, calculation


@pytest.fixture(scope="module")
def mock_qeapp_workchain():
    qe_workchain = mock_workchain("QeAppWorkChain")
    relax_workchain = mock_workchain("PwRelaxWorkChain")
    base_workchain = mock_workchain("PwBaseWorkChain")
    relax_calcjob = mock_calcjob("PwCalculation")
    parameters = orm.Dict(dict={"CONTROL": {"calculation": "relax"}})
    parameters.store()
    relax_calcjob.base.links.add_incoming(
        parameters,
        link_type=LinkType.INPUT_CALC,
        link_label="parameters",
    )
    relax_calcjob.base.links.add_incoming(
        base_workchain,
        link_type=LinkType.CALL_CALC,
        link_label="iteration_01",
    )
    base_workchain.base.links.add_incoming(
        parameters,
        link_type=LinkType.INPUT_WORK,
        link_label="pw__parameters",
    )
    base_workchain.base.links.add_incoming(
        relax_workchain,
        link_type=LinkType.CALL_WORK,
        link_label="iteration_01",
    )
    relax_workchain.base.links.add_incoming(
        qe_workchain,
        link_type=LinkType.CALL_WORK,
        link_label="relax",
    )
    base_workchain.set_metadata_inputs(
        {
            "pw": {"metadata": {"options": {}}},
            "metadata": {},
        }
    )
    relax_workchain.set_metadata_inputs(
        {
            "base": {
                "pw": {"metadata": {"options": {}}},
                "metadata": {},
            },
            "metadata": {},
        }
    )
    qe_workchain.set_metadata_inputs(
        {
            "relax": {
                "base": {
                    "pw": {"metadata": {"options": {}}},
                },
            },
        }
    )
    qe_workchain.store()
    relax_workchain.store()
    base_workchain.store()
    relax_calcjob.store()
    yield qe_workchain


class TreeTestingMixin:
    tree: SimplifiedProcessTree

    @property
    def relax_workchain_node(self) -> WorkChainTreeNode:
        return self.tree.trunk.branches[0]  # type: ignore

    @property
    def base_workchain_node(self) -> WorkChainTreeNode:
        return self.relax_workchain_node.branches[0]  # type: ignore

    @property
    def calcjob_node(self) -> CalcJobTreeNode:
        return self.base_workchain_node.branches[0]  # type: ignore


class TestSimplifiedProcessTree(TreeTestingMixin):
    model: SimplifiedProcessTreeModel
    node: orm.WorkChainNode

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, request):
        model = SimplifiedProcessTreeModel()
        request.cls.model = model
        request.cls.tree = SimplifiedProcessTree(model=model)

    def test_initialization(self):
        loading_widget = self.tree.children[0]  # type: ignore
        assert loading_widget is self.tree.loading_message
        assert (
            loading_widget.message.value
            == self.tree.loading_message.message.value
            == "Loading process tree"
        )
        assert "simplified-process-tree" in self.tree._dom_classes
        assert not self.tree.rendered

    def test_render(self):
        self.tree.render()
        assert self.tree.children[0].value == "Process node not yet available"
        assert not self.tree.rendered  # no process node yet

    def test_update_on_process_change(self, mock_qeapp_workchain):
        self.model.process_uuid = mock_qeapp_workchain.uuid
        assert self.tree.rendered
        human_label = TITLE_MAPPING["QeAppWorkChain"]
        assert self.tree.trunk.label.value == human_label
        assert not self.tree.trunk.collapsed
        assert len(self.tree.trunk.branches) == 1

    def test_relax_workchain_manager_node(self):
        assert isinstance(self.relax_workchain_node, WorkChainTreeNode)
        assert self.relax_workchain_node.level == 1
        assert self.relax_workchain_node.emoji.value == "✅"
        assert self.relax_workchain_node.state.value == "finished"
        human_label = TITLE_MAPPING["PwRelaxWorkChain"]
        assert self.relax_workchain_node.label.value == human_label
        assert isinstance(self.relax_workchain_node.label, ipw.HTML)
        assert self.relax_workchain_node.collapsed
        assert len(self.relax_workchain_node.branches) == 0

    def test_expand(self):
        self.relax_workchain_node.expand()
        assert "open" in self.relax_workchain_node.branches._dom_classes
        assert self.relax_workchain_node.toggle.icon == "minus"
        assert len(self.relax_workchain_node.branches) == 1

    def test_collapse(self):
        self.relax_workchain_node.collapse()
        assert "open" not in self.relax_workchain_node.branches._dom_classes
        assert self.relax_workchain_node.toggle.icon == "plus"

    def test_expand_recursive(self):
        self.tree.trunk.expand(recursive=True)
        assert all(
            not branch.collapsed
            for branch in (
                self.tree.trunk,
                self.relax_workchain_node,
                self.base_workchain_node,
            )
        )

    def test_workchain_node(self):
        assert isinstance(self.base_workchain_node, WorkChainTreeNode)
        assert self.base_workchain_node.level == 2
        assert self.base_workchain_node.emoji.value == "✅"
        assert self.base_workchain_node.state.value == "finished"
        human_label = TITLE_MAPPING["PwBaseWorkChain"]["relax"]
        assert self.base_workchain_node.label.value == human_label
        assert isinstance(self.base_workchain_node.label, ipw.HTML)
        assert len(self.base_workchain_node.branches) == 1

    def test_collapse_recursive(self):
        self.relax_workchain_node.collapse(recursive=True)
        assert self.relax_workchain_node.collapsed
        assert self.base_workchain_node.collapsed

    def test_collapse_all_button(self):
        self.tree.trunk.expand(recursive=True)
        self.tree.collapse_button.click()
        assert all(
            branch.collapsed
            for branch in (
                self.tree.trunk,
                self.relax_workchain_node,
                self.base_workchain_node,
            )
        )

    def test_calcjob_node(self):
        assert isinstance(self.calcjob_node, CalcJobTreeNode)
        assert self.calcjob_node.level == 3
        assert self.calcjob_node.emoji.value == "✅"
        assert self.calcjob_node.state.value == "finished"
        human_label = TITLE_MAPPING["PwCalculation"]["relax"]
        assert isinstance(self.calcjob_node.label, ipw.Button)
        assert self.calcjob_node.label.description == human_label

    @pytest.mark.parametrize("without_flag", [True, False])
    def test_tree_monitor_updates(self, without_flag):
        self.tree.trunk.collapse()
        self.tree.trunk.clear()
        assert not self.tree.trunk.branches

        def update_monitor():
            for _ in range(10):
                self.model.monitor_counter += 1

        if without_flag:  # skip `_adding_branches` flag check
            add_branches = self.tree.trunk._add_branches  # store original method

            def delayed_add_branches_recursive(
                node: orm.ProcessNode | None = None,
            ):
                node = node or self.tree.trunk.process
                for child in sorted(node.called, key=lambda child: child.ctime):
                    if isinstance(child, orm.CalcFunctionNode):
                        continue
                    if child.pk in self.tree.trunk.pks:
                        continue
                    if child.process_label in (
                        "BandsWorkChain",
                        "ProjwfcBaseWorkChain",
                    ):
                        delayed_add_branches_recursive(child)
                    else:
                        time.sleep(0.5)
                        self.tree.trunk._add_branch(child)

            self.tree.trunk._add_branches = delayed_add_branches_recursive

        # Start monitoring thread to update the tree, which will also
        # attempt to add branches to expanded parents
        monitor_thread = threading.Thread(target=update_monitor)
        monitor_thread.start()

        # Expand the trunk in the main thread
        self.tree.trunk.toggle.click()

        # Wait for the background thread to finish before asserting
        monitor_thread.join()

        if without_flag:
            # Without the `_adding_branches` flag, the monitor thread will
            # attempt to add branches at the same time as the main thread,
            # which will result in duplicate branches
            assert len(self.tree.trunk.branches) > 1
            self.tree.trunk._add_branches = add_branches  # restore original method
        else:
            # With the flag in place, the monitor thread bails from adding
            # branches when the parent branch is presently doing so, leading
            # to only the adding of the single intended branch
            assert len(self.tree.trunk.branches) == 1


class TestWorkChainStatusPanel(TreeTestingMixin):
    model: WorkChainStatusModel
    panel: WorkChainStatusPanel

    @property
    def tree(self) -> SimplifiedProcessTree:
        return self.panel.simplified_process_tree

    @pytest.fixture(scope="class", autouse=True)
    def setup(self, request, mock_qeapp_workchain):
        model = WorkChainStatusModel()
        model.process_uuid = mock_qeapp_workchain.uuid
        request.cls.model = model
        request.cls.panel = WorkChainStatusPanel(model=model)

    def test_render(self):
        self.panel.render()
        assert self.panel.children[0] is self.panel.accordion
        assert self.panel.accordion.selected_index == 1
        assert self.panel.simplified_process_tree.rendered

    def test_calcjob_node_link(self):
        # TODO uncomment if and when the comments below are resolved; discard otherwise
        # from aiidalab_qe.common.node_view import CalcJobNodeViewerWidget

        trunk = self.tree.trunk
        trunk.expand(recursive=True)
        self.calcjob_node.label.click()
        assert self.panel.process_tree.value == self.calcjob_node.process_uuid
        assert self.panel.process_tree._tree.nodes == (self.calcjob_node.process,)
        # TODO understand why the following does not trigger automatically as in the app
        # TODO understand why the following triggers a thread
        # self.panel.process_tree.set_trait("selected_nodes", [self.calcjob_node.node])
        # assert isinstance(self.panel.node_view, CalcJobNodeViewerWidget)
        # assert self.panel.node_view_container.children[0] is self.node_view  # type: ignore

    def test_to_advanced_view_button(self):
        assert self.panel.accordion.selected_index == 1
        self.panel.to_advanced_view_button.click()
        assert self.panel.accordion.selected_index == 2

    def test_paused_title_and_navigation(self, mock_qeapp_workchain):
        calculation = next(
            node
            for node in mock_qeapp_workchain.called_descendants
            if isinstance(node, orm.CalcJobNode)
        )
        self.panel.paused_processes_model.paused_count = 2
        assert self.panel.accordion.titles[0] == "Paused Processes (2)"

        self.panel._on_paused_inspect(calculation.uuid)

        assert self.panel.accordion.selected_index == 2
        assert self.panel.process_tree.value == calculation.uuid


def test_paused_model_detects_paused_root_and_descendants():
    root, child, calculation = create_process_graph()
    root.pause()
    root.set_process_status("Root paused")
    calculation.pause()
    calculation.set_process_status("Calculation paused")

    model = PausedProcessesModel()
    model.process_uuid = root.uuid
    model.update()

    assert model.paused_count == 2
    assert {node.uuid for node in model.nodes} == {root.uuid, calculation.uuid}
    assert child.uuid not in {node.uuid for node in model.nodes}


def test_paused_model_updates_when_monitor_counter_changes():
    root, _, calculation = create_process_graph()
    model = PausedProcessesModel()
    model.process_uuid = root.uuid

    model.monitor_counter += 1
    assert model.paused_count == 0

    calculation.pause()
    model.monitor_counter += 1
    assert model.paused_count == 1
    assert model.nodes[0].uuid == calculation.uuid


def test_paused_model_play_clears_previous_error(monkeypatch):
    root, _, _ = create_process_graph()
    root.pause()
    model = PausedProcessesModel()
    model.process_uuid = root.uuid
    model.error_message = "previous error"
    play_processes = Mock()
    monkeypatch.setattr(control, "play_processes", play_processes)

    model.play(root.uuid)

    play_processes.assert_called_once()
    assert play_processes.call_args.args[0][0].uuid == root.uuid
    assert model.error_message == ""


def test_paused_model_play_reports_errors(monkeypatch):
    root, _, _ = create_process_graph()
    model = PausedProcessesModel()
    model.process_uuid = root.uuid
    monkeypatch.setattr(
        control,
        "play_processes",
        Mock(side_effect=RuntimeError("daemon unavailable")),
    )

    model.play(root.uuid)

    assert model.error_message == "daemon unavailable"


class TestPausedProcessesTable:
    def test_renders_rows_and_actions(self):
        root, _, calculation = create_process_graph()
        calculation.pause()
        calculation.set_process_status("SCF paused")

        model = PausedProcessesModel()
        model.process_uuid = root.uuid
        model.update()
        on_inspect = Mock()
        panel = PausedProcessesTable(model=model, on_inspect=on_inspect)

        assert len(panel.table.children) == 2
        assert panel.table.children[0].children[0].value == "<b>PK</b>"

        row = panel.table.children[1]
        assert len(row.children) == 5
        assert row.children[0].value == f"<b>{calculation.pk}</b>"
        assert row.children[2].value == "SCF paused"

        row.children[3].click()
        on_inspect.assert_called_once_with(calculation.uuid)

    def test_play_button_delegates_to_model(self, monkeypatch):
        root, _, calculation = create_process_graph()
        calculation.pause()
        model = PausedProcessesModel()
        model.process_uuid = root.uuid
        model.update()
        play = Mock()
        monkeypatch.setattr(model, "play", play)
        panel = PausedProcessesTable(model=model)

        panel.table.children[1].children[4].click()

        play.assert_called_once_with(calculation.uuid)

    def test_displays_error_message(self):
        model = PausedProcessesModel()
        panel = PausedProcessesTable(model=model)

        model.error_message = "daemon unavailable"

        assert "alert-danger" in panel.alert.value
        assert "daemon unavailable" in panel.alert.value

        model.error_message = ""
        assert panel.alert.value == ""

    def test_daemon_warning_disables_play(self):
        root, _, calculation = create_process_graph()
        model = PausedProcessesModel()
        model.process_uuid = root.uuid
        model.update()
        model.daemon_is_running = False
        panel = PausedProcessesTable(model=model)

        assert panel.daemon_warning.value == ""

        calculation.pause()
        model.monitor_counter += 1

        assert "daemon is not running" in panel.daemon_warning.value
        assert panel.table.children[1].children[4].disabled

        model.daemon_is_running = True

        assert panel.daemon_warning.value == ""
        assert not panel.table.children[1].children[4].disabled
