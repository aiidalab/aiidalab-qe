def test_result_step(app_to_submit, generate_qeapp_workchain):
    """Test the result step is properly updated when the process
    is running."""

    step = app_to_submit.results_step
    step.process = generate_qeapp_workchain().node.uuid
    assert step.state == step.State.ACTIVE


def test_workchainview(generate_qeapp_workchain):
    """Test the result tabs are properly updated"""
    from aiidalab_qe.app.result.workchain_viewer import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wcv = WorkChainViewer(wkchain.node)
    assert len(wcv.result_tabs.children) == 5
    assert wcv.result_tabs._titles["0"] == "Workflow Summary"
    assert wcv.result_tabs._titles["1"] == "Final Geometry"
    assert wcv.result_tabs._titles["2"] == "Electronic Structure"


def test_electronic_structure(generate_qeapp_workchain):
    """test the report can be properly generated from the builder without errors"""
    from aiida import engine

    from aiidalab_qe.app.result.workchain_viewer import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wkchain.node.set_exit_status(0)
    wkchain.node.set_process_state(engine.ProcessState.FINISHED)
    wcv = WorkChainViewer(wkchain.node)
    assert (
        "DOS grouped by:"
        == wcv.result_tabs.children[2].children[0].children[0].children[0].value
    )


def test_summary_report(data_regression, generate_qeapp_workchain):
    """Test the summary report can be properly generated."""
    from aiidalab_qe.app.result.summary_viewer import SummaryView

    wkchain = generate_qeapp_workchain()
    viewer = SummaryView(wkchain.node)
    report = viewer.report
    # regression test
    data_regression.check(report)


def test_summary_report_advanced_settings(data_regression, generate_qeapp_workchain):
    """Test advanced settings are properly reported"""
    from aiidalab_qe.app.result.summary_viewer import SummaryView

    wkchain = generate_qeapp_workchain(
        spin_type="collinear", initial_magnetic_moments=0.1
    )
    viewer = SummaryView(wkchain.node)
    report = viewer.report
    assert report["initial_magnetic_moments"]["Si"] == 0.1


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
