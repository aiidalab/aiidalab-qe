def test_result_step(app_to_submit, generate_qeapp_workchain):
    """Test the result step is properly updated when the process
    is running."""

    step = app_to_submit.steps.steps[3][1]
    step.process = generate_qeapp_workchain().node.uuid
    assert step.state == step.State.ACTIVE


def test_workchainview(generate_qeapp_workchain):
    from aiidalab_qe.app.result.node_view import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wcv = WorkChainViewer(wkchain.node)
    assert len(wcv.result_tabs.children) == 5
    assert wcv.result_tabs._titles["0"] == "Workflow Summary"
    assert wcv.result_tabs._titles["1"] == "Final Geometry"
    assert wcv.result_tabs._titles["2"] == "Electronic Structure"


def test_summary_view(generate_qeapp_workchain):
    """test the report can be properly generated from the builder without errors"""
    from bs4 import BeautifulSoup

    from aiidalab_qe.app.result.node_view import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wcv = WorkChainViewer(wkchain.node)
    report_html = wcv.result_tabs.children[0].children[0].children[0].value
    # report_html = generate_report_html(wcv.node)
    parsed = BeautifulSoup(report_html, "html.parser")
    print(parsed)
    # find the td with the text "Initial Magnetic Moments"
    parameters = {
        "Energy cutoff (wave functions)": "30.0 Ry",
        "Total Charge": "0.0",
        "Initial Magnetic Moments": "None",
    }
    for key, value in parameters.items():
        td = parsed.find("td", text=key).find_next_sibling("td")
        assert td.text == value


def test_electronic_structure(generate_qeapp_workchain):
    """test the report can be properly generated from the builder without errors"""
    from aiidalab_qe.app.result.node_view import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wcv = WorkChainViewer(wkchain.node)
    print(wcv.result_tabs.children[2].children[0].children[0].children[0].value)
    assert (
        "DOS grouped by:"
        == wcv.result_tabs.children[2].children[0].children[0].children[0].value
    )


def test_workchain_outputs(generate_qeapp_workchain):
    """Test the workchain outputs widget can be properly generated
    without errors"""
    from aiida import engine

    from aiidalab_qe.app.result.node_view import WorkChainOutputs

    wkchain = generate_qeapp_workchain()
    qeapp_node = wkchain.node
    qeapp_node.set_exit_status(0)
    qeapp_node.set_process_state(engine.ProcessState.FINISHED)
    widget = WorkChainOutputs(qeapp_node)
    widget._download_archive(1)
    # TODO: check the content of the archive
