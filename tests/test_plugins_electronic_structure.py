def test_electronic_structure(generate_qeapp_workchain):
    """Test the electronic structure tab."""
    from aiida import engine

    from aiidalab_qe.app.result.workchain_viewer import WorkChainViewer

    wkchain = generate_qeapp_workchain()
    wkchain.node.set_exit_status(0)
    wkchain.node.set_process_state(engine.ProcessState.FINISHED)
    wcv = WorkChainViewer(wkchain.node)
    # find the tab with the identifier "electronic_structure"
    # the built-in summary and structure tabs is not a plugin panel,
    # thus don't have identifiers
    tab = [
        tab
        for tab in wcv.result_tabs.children
        if getattr(tab, "identifier", "") == "electronic_structure"
    ][0]
    # It should have two children: settings and the _bands_plot_view
    assert len(tab.children) == 2
