def test_workchainview(generate_qeapp_workchain):
    from aiidalab_qe.app.result.node_view import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wcv = WorkChainViewer(wkchain.node)
    assert len(wcv.result_tabs.children) == 5
    assert wcv.result_tabs._titles["0"] == "Workflow Summary"
    assert wcv.result_tabs._titles["1"] == "Final Geometry"
    assert wcv.result_tabs._titles["2"] == "Electronic Structure"
