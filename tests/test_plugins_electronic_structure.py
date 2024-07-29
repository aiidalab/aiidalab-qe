def test_electronic_structure(generate_qeapp_workchain):
    """Test the electronic structure tab."""
    import plotly.graph_objects as go
    from aiidalab_qe.app.result.workchain_viewer import WorkChainViewer
    from aiidalab_qe.common.bandpdoswidget import BandPdosWidget
    from aiidalab_qe.plugins.electronic_structure.result import Result

    from aiida import engine

    wkchain = generate_qeapp_workchain()
    wkchain.node.set_exit_status(0)
    wkchain.node.set_process_state(engine.ProcessState.FINISHED)
    wcv = WorkChainViewer(wkchain.node)
    # find the tab with the identifier "electronic_structure"
    # the built-in summary and structure tabs is not a plugin panel,
    # thus don't have identifiers
    tab = next(
        tab
        for tab in wcv.result_tabs.children
        if getattr(tab, "identifier", "") == "electronic_structure"
    )
    # It should have one children: the _bands_plot_view
    assert len(tab.children) == 1

    result = Result(node=wkchain.node)
    result._update_view()

    assert isinstance(result.children[0], BandPdosWidget)
    assert isinstance(result.children[0].bandsplot_widget, go.FigureWidget)

    # Check if data is correct
    assert result.children[0].bands_data is not None
    assert result.children[0].bands_data["pathlabels"] is not None
    assert result.children[0].pdos_data is not None

    # Check Bands axis
    assert result.children[0].bandsplot_widget.layout.xaxis.title.text == "k-points"
    assert (
        result.children[0].bandsplot_widget.layout.xaxis2.title.text
        == "Density of states"
    )
    assert (
        result.children[0].bandsplot_widget.layout.yaxis.title.text
        == "Electronic Bands (eV)"
    )
    assert isinstance(
        result.children[0].bandsplot_widget.layout.xaxis.rangeslider,
        go.layout.xaxis.Rangeslider,
    )
    assert result.children[0].bands_data["pathlabels"][0] == list(
        result.children[0].bandsplot_widget.layout.xaxis.ticktext
    )
