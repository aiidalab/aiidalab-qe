import pytest


@pytest.mark.usefixtures("sssp")
def test_result_step(app_to_submit, generate_qeapp_workchain):
    """Test the result step is properly updated when the process
    is running."""

    step = app_to_submit.results_step
    step.process = generate_qeapp_workchain().node.uuid
    assert step.state == step.State.ACTIVE


@pytest.mark.usefixtures("sssp")
def test_kill_and_clean_buttons(app_to_submit, generate_qeapp_workchain):
    """Test the kill and clean_scratch button are properly displayed when the process
    is in different states."""

    step = app_to_submit.results_step
    step.process = generate_qeapp_workchain().node.uuid
    step._update_state()
    step._update_kill_button_layout()
    step._update_clean_scratch_button_layout()
    assert step.kill_button.layout.display == "block"
    assert step.clean_scratch_button.layout.display == "none"


@pytest.mark.usefixtures("sssp")
def test_workchainview(generate_qeapp_workchain):
    """Test the result tabs are properly updated"""
    import time

    from aiidalab_qe.app.result.workchain_viewer import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wcv = WorkChainViewer(wkchain.node)
    # wait for the tabs to be updated by the process monitor
    time.sleep(3)
    assert len(wcv.result_tabs.children) == 5
    assert wcv.result_tabs._titles["0"] == "Workflow Summary"
    assert wcv.result_tabs._titles["1"] == "Final Geometry"


@pytest.mark.usefixtures("sssp")
def test_summary_report(data_regression, generate_qeapp_workchain):
    """Test the summary report can be properly generated."""
    from aiidalab_qe.app.result.summary_viewer import SummaryView

    wkchain = generate_qeapp_workchain()
    viewer = SummaryView(wkchain.node)
    report = viewer.report
    # regression test
    data_regression.check(report)


@pytest.mark.usefixtures("sssp")
def test_summary_report_advanced_settings(data_regression, generate_qeapp_workchain):
    """Test advanced settings are properly reported"""
    from aiidalab_qe.app.result.summary_viewer import SummaryView

    wkchain = generate_qeapp_workchain(
        spin_type="collinear", electronic_type="metal", initial_magnetic_moments=0.1
    )
    viewer = SummaryView(wkchain.node)
    report = viewer.report
    assert report["initial_magnetic_moments"]["Si"] == 0.1


@pytest.mark.usefixtures("sssp")
def test_summary_view(generate_qeapp_workchain):
    """Test the report html can be properly generated."""
    from bs4 import BeautifulSoup

    from aiidalab_qe.app.result.summary_viewer import SummaryView

    wkchain = generate_qeapp_workchain()
    viewer = SummaryView(wkchain.node)
    report_html = viewer.report_html
    # report_html = generate_report_html(wcv.node)
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
