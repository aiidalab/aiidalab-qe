def test_electronic_structure(generate_qeapp_workchain):
    """test the report can be properly generated from the builder without errors"""
    from aiida import engine

    from aiidalab_qe.app.result.workchain_viewer import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wkchain.node.set_exit_status(0)
    wkchain.node.set_process_state(engine.ProcessState.FINISHED)
    wcv = WorkChainViewer(wkchain.node)
    # find the tab with the identifier "electronic_structure"
    tab = [
        tab
        for tab in wcv.result_tabs.children
        if tab.identifier == "electronic_structure"
    ][0]
    # check the content of the tab
    assert "DOS grouped by:" == tab.children[0].children[0].value
