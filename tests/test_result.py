import typing as t
from threading import Event
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from bs4 import BeautifulSoup

from aiida.engine.processes import control
from aiidalab_qe.app.result import ResultsStep, ResultsStepModel
from aiidalab_qe.app.result.components.status.paused import PausedProcess
from aiidalab_qe.app.result.components.summary import WorkflowSummaryModel
from aiidalab_qe.app.result.components.viewer import (
    WorkflowResultsViewer,
    WorkflowResultsViewerModel,
)
from aiidalab_qe.app.result.components.viewer.structure import (
    StructureResultsModel,
    StructureResultsPanel,
)
from aiidalab_qe.app.wizard import QeWizard
from aiidalab_qe.common import widgets as common_widgets
from aiidalab_qe.common.wizard import State


def test_calcjob_output_follower_threads_are_daemon(monkeypatch):
    started = Event()
    monkeypatch.setattr(
        common_widgets, "load_node", lambda _: SimpleNamespace(is_sealed=False)
    )
    follower = common_widgets.CalcJobOutputFollower()

    def fetch_output(_):
        started.set()
        return []

    monkeypatch.setattr(follower, "_fetch_output", fetch_output)
    follower.calcjob_uuid = "calcjob-uuid"
    try:
        assert started.wait(timeout=2)
        assert follower._follow_output_thread.daemon
        assert follower._push_thread.is_alive() and follower._push_thread.daemon
        assert follower._pull_thread.is_alive() and follower._pull_thread.daemon
    finally:
        follower._stop_follow_output.set()
        for thread in (
            follower._follow_output_thread,
            follower._push_thread,
            follower._pull_thread,
        ):
            if thread is not None:
                thread.join(timeout=2)
    assert not follower._push_thread.is_alive()
    assert not follower._pull_thread.is_alive()


def test_result_step(app_to_submit, generate_qeapp_workchain):
    """Test the result step is properly updated when the process
    is running."""
    app: QeWizard = app_to_submit
    step: ResultsStep = app.results_step
    model: ResultsStepModel = app.results_model
    model.process_uuid = generate_qeapp_workchain().node.uuid
    assert model.state == State.ACTIVE
    step.render()
    assert step.toggle_controls.value == "Status"
    results_viewer_model = t.cast(
        WorkflowResultsViewerModel, model.get_model("results")
    )
    # All jobs are completed, so there should be no process status notifications
    for _, results_model in results_viewer_model.get_models():
        assert results_model.process_status_notification == ""


def test_kill_and_clean_buttons(app_to_submit, generate_qeapp_workchain):
    """Test the kill and clean_scratch button are properly displayed when the process
    is in different states."""
    app: QeWizard = app_to_submit
    model: ResultsStepModel = app.results_model
    step: ResultsStep = app.results_step
    step.render()
    model.process_uuid = generate_qeapp_workchain().node.uuid
    assert step.kill_button.layout.display == "block"
    assert step.clean_scratch_button.layout.display == "none"


def test_paused_processes_warning(app_to_submit):
    app: QeWizard = app_to_submit
    step: ResultsStep = app.results_step
    step.render()

    assert step.paused_processes_warning.value == ""

    step.status_panel.paused_processes_model.paused_count = 2

    assert "2 paused processes" in step.paused_processes_warning.value
    assert "Paused processes" in step.paused_processes_warning.value


def test_kill_workflow_resets_paused_processes(
    app_to_submit,
    generate_qeapp_workchain,
    monkeypatch,
):
    app: QeWizard = app_to_submit
    step: ResultsStep = app.results_step
    model: ResultsStepModel = app.results_model
    model.process_uuid = generate_qeapp_workchain().node.uuid
    step.render()
    paused_model = step.status_panel.paused_processes_model
    model.daemon_is_running = True
    model.daemon_status_known = True
    paused_model.daemon_is_running = True
    paused_model.daemon_status_known = True
    paused_model.paused_processes = (
        PausedProcess(
            uuid=model.process.uuid,
            pk=model.process.pk,
            label=model.process.label,
            status="Paused",
        ),
    )
    paused_model.paused_count = 1
    monkeypatch.setattr(control, "kill_processes", Mock())

    model.kill_pending = True
    model._kill_deadline = 0
    paused_model.paused_processes = ()
    model.process_pending_kill()
    paused_model.reset()

    assert paused_model.paused_processes == ()
    assert paused_model.paused_count == 0


def test_pending_kill_failure_keeps_monitoring_and_allows_retry(
    app_to_submit,
    generate_qeapp_workchain,
    monkeypatch,
):
    app: QeWizard = app_to_submit
    step: ResultsStep = app.results_step
    model: ResultsStepModel = app.results_model
    model.process_uuid = generate_qeapp_workchain().node.uuid
    step.render()
    model.daemon_is_running = True
    model.daemon_status_known = True
    model.update_daemon_status = Mock()
    model.kill_pending = True
    model._kill_deadline = 0
    kill_processes = Mock(side_effect=[RuntimeError("kill RPC failed"), None])
    monkeypatch.setattr(control, "kill_processes", kill_processes)
    initial_monitor_counter = model.monitor_counter

    step._update_status()

    assert model.monitor_counter > initial_monitor_counter
    assert model.kill_pending is False
    assert "kill RPC failed" in step.status_panel.paused_processes_model.error_message

    model.kill_process()

    assert kill_processes.call_count == 2


def test_workchainview(generate_qeapp_workchain):
    """Test the result tabs are properly updated"""
    workchain = generate_qeapp_workchain()
    workchain.node.seal()
    model = WorkflowResultsViewerModel()
    viewer = WorkflowResultsViewer(model=model)
    model.process_uuid = workchain.node.uuid
    viewer.render()
    assert len(viewer.tabs.children) == 2
    assert viewer.tabs.titles[0] == "Structure"


def test_summary_report(data_regression, generate_qeapp_workchain):
    """Test the summary report can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkflowSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    # Discard variable parameters
    for key in (
        "pk",
        "uuid",
        "creation_time",
        "modification_time",
    ):
        report_parameters["workflow_properties"].pop(key)
    for key in (
        "structure_pk",
        "structure_uuid",
    ):
        report_parameters["initial_structure_properties"].pop(key)
    data_regression.check(report_parameters)


def test_summary_report_advanced_settings(data_regression, generate_qeapp_workchain):
    """Test advanced settings are properly reported"""
    workchain = generate_qeapp_workchain(
        spin_type="collinear", electronic_type="metal", initial_magnetic_moments=0.1
    )
    model = WorkflowSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    moments = report_parameters["advanced_settings"]["initial_magnetic_moments"]
    assert moments["Si"] == 0.1


@pytest.mark.parametrize(
    ("pbc", "symmetry_key"),
    [
        [(False, False, False), "point_group"],  # 0D
        [(True, False, False), "space_group"],  # 1D
        [(True, True, False), "space_group"],  # 2D
        [(True, True, True), "space_group"],  # 3D
    ],
)
def test_summary_report_symmetry_group(
    generate_qeapp_workchain,
    generate_structure_data,
    pbc,
    symmetry_key,
):
    """Test summary report includes correct symmetry group for all system dimension."""

    system = generate_structure_data("silicon", pbc=pbc)
    workchain = generate_qeapp_workchain(
        structure=system,
        run_bands=False,
        relax_type="none",
    )
    model = WorkflowSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    assert symmetry_key in report_parameters["initial_structure_properties"]


def test_summary_view(generate_qeapp_workchain):
    """Test the report html can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkflowSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_html = model.generate_report_html()
    parsed = BeautifulSoup(report_html, "html.parser")
    parameters = {
        "Energy cutoff (wave functions)": "30.0 Ry",
        "Total charge": "0.0",
    }
    for key, value in parameters.items():
        key_td = parsed.find("td", string=lambda tag, key=key: tag and key in tag.text)
        value_td = key_td.find_next_sibling("td")
        assert value in value_td.text


def test_structure_results_panel(generate_qeapp_workchain):
    """Test the structure results panel can be properly generated."""

    model = StructureResultsModel()
    panel = StructureResultsPanel(model=model)

    def test_table_data(model: StructureResultsModel):
        rows = model.table_data[1:]  # skip table header # type: ignore
        for i, row in enumerate(rows):
            position = model.structure.sites[i].position
            x, y, z = (f"{coordinate:.2f}" for coordinate in position)  # type: ignore
            assert row == [i + 1, "Si", 0, x, y, z]  # type: ignore

    assert model.title == "Structure"

    wc = generate_qeapp_workchain(relax_type="none")
    model.process_uuid = wc.node.uuid
    assert "Si<sub>2</sub>" in model.header
    assert "Initial" in model.sub_header
    assert "properties" in model.source  # inputs
    assert model.structure.pk == model.inputs.structure.pk
    assert str(model.inputs.structure.pk) in model.info
    test_table_data(model)

    panel.render()
    assert panel.view_toggle_button.layout.display == "none"

    wc = generate_qeapp_workchain(relax_type="positions_cell")
    model.process_uuid = wc.node.uuid
    assert "Initial" in model.sub_header
    assert panel.view_toggle_button.layout.display == "block"
    assert panel.view_toggle_button.description == "View relaxed"
    panel.view_toggle_button.click()
    assert panel.view_toggle_button.description == "View initial"
    assert "Relaxed" in model.sub_header
    assert "properties" not in model.source  # outputs
    assert model.structure.pk == model.outputs.structure.pk
    assert str(model.outputs.structure.pk) in model.info
    test_table_data(model)
