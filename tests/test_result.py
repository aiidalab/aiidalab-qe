import time

from bs4 import BeautifulSoup

from aiidalab_qe.app.main import App
from aiidalab_qe.app.result.summary import WorkChainSummaryModel
from aiidalab_qe.app.result.viewer import WorkChainViewer, WorkChainViewerModel


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
    model = WorkChainViewerModel()
    viewer = WorkChainViewer(workchain.node, model=model)
    time.sleep(3)
    assert len(viewer.tabs.children) == 5
    assert viewer.tabs._titles["0"] == "Workflow Summary"  # type: ignore
    assert viewer.tabs._titles["1"] == "Final Geometry"  # type: ignore


def test_summary_report(data_regression, generate_qeapp_workchain):
    """Test the summary report can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkChainSummaryModel()
    model.process_node = workchain.node
    report_parameters = model._generate_report_parameters()
    data_regression.check(report_parameters)


def test_summary_report_advanced_settings(data_regression, generate_qeapp_workchain):
    """Test advanced settings are properly reported"""
    workchain = generate_qeapp_workchain(
        spin_type="collinear", electronic_type="metal", initial_magnetic_moments=0.1
    )
    model = WorkChainSummaryModel()
    model.process_node = workchain.node
    report_parameters = model._generate_report_parameters()
    assert report_parameters["initial_magnetic_moments"]["Si"] == 0.1


def test_summary_view(generate_qeapp_workchain):
    """Test the report html can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkChainSummaryModel()
    model.process_node = workchain.node
    report_html = model.generate_report_html()
    parsed = BeautifulSoup(report_html, "html.parser")
    # find the td with the text "Initial Magnetic Moments"
    parameters = {
        "Energy cutoff (wave functions)": "30.0 Ry",
        "Total Charge": "0.0",
        "Initial Magnetic Moments": "",
    }
    for key, value in parameters.items():
        td = parsed.find("td", text=key).find_next_sibling("td")
        assert td.text == value
