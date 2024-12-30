from bs4 import BeautifulSoup

from aiidalab_qe.app.main import App
from aiidalab_qe.app.result.components.summary import WorkChainSummaryModel
from aiidalab_qe.app.result.components.viewer import (
    WorkChainResultsViewer,
    WorkChainResultsViewerModel,
)
from aiidalab_qe.app.result.components.viewer.structure import StructureResultsModel
from aiidalab_qe.app.result.components.viewer.structure.structure import (
    StructureResultsPanel,
)


def test_result_step(app_to_submit, generate_qeapp_workchain):
    """Test the result step is properly updated when the process
    is running."""
    app: App = app_to_submit
    step = app.results_step
    app.results_model.process_uuid = generate_qeapp_workchain().node.uuid
    assert step.state == step.State.ACTIVE


def test_kill_and_clean_buttons(app_to_submit, generate_qeapp_workchain):
    """Test the kill and clean_scratch button are properly displayed when the process
    is in different states."""
    step = app_to_submit.results_step
    step.render()
    model = app_to_submit.results_model
    model.process_uuid = generate_qeapp_workchain().node.uuid
    assert step.kill_button.layout.display == "block"
    assert step.clean_scratch_button.layout.display == "none"


def test_workchainview(generate_qeapp_workchain):
    """Test the result tabs are properly updated"""
    workchain = generate_qeapp_workchain()
    workchain.node.seal()
    model = WorkChainResultsViewerModel()
    viewer = WorkChainResultsViewer(model=model)
    model.process_uuid = workchain.node.uuid
    viewer.render()
    assert len(viewer.tabs.children) == 4
    assert viewer.tabs._titles["0"] == "Relaxed structure"  # type: ignore


def test_summary_report(data_regression, generate_qeapp_workchain):
    """Test the summary report can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkChainSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    # Discard variable parameters
    for key in (
        "pk",
        "uuid",
        "creation_time",
        "creation_time_relative",
        "modification_time",
        "modification_time_relative",
    ):
        report_parameters.pop(key)
    data_regression.check(report_parameters)


def test_summary_report_advanced_settings(data_regression, generate_qeapp_workchain):
    """Test advanced settings are properly reported"""
    workchain = generate_qeapp_workchain(
        spin_type="collinear", electronic_type="metal", initial_magnetic_moments=0.1
    )
    model = WorkChainSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    assert report_parameters["initial_magnetic_moments"]["Si"] == 0.1


def test_summary_view(generate_qeapp_workchain):
    """Test the report html can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkChainSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_html = model.generate_report_html()
    parsed = BeautifulSoup(report_html, "html.parser")
    parameters = {
        "Energy cutoff (wave functions)": "30.0 Ry",
        "Total charge": "0.0",
    }
    for key, value in parameters.items():
        li = parsed.find_all(lambda tag, key=key: tag.name == "li" and key in tag.text)
        assert value in li[0].text


def test_structure_results_panel(generate_qeapp_workchain):
    """Test the structure results panel can be properly generated."""
    model = StructureResultsModel()
    _ = StructureResultsPanel(model=model)

    wc = generate_qeapp_workchain(relax_type="none")
    model.process_uuid = wc.node.uuid
    assert model.title == "Initial structure"
    assert "properties" in model.source  # source should be inputs

    wc = generate_qeapp_workchain(relax_type="positions_cell")
    model.process_uuid = wc.node.uuid
    assert model.title == "Relaxed structure"
    assert "properties" not in model.source  # source should be outputs
