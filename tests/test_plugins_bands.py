def test_result(generate_qeapp_workchain):
    from widget_bandsplot import BandsPlotWidget

    from aiidalab_qe.plugins.bands.result import Result, export_bands_data

    wkchain = generate_qeapp_workchain()
    #
    data = export_bands_data(wkchain.node.outputs.bands)
    assert data is not None
    # generate structure for scf calculation
    result = Result(wkchain.node)
    assert result.identifier == "bands"
    result._update_view()
    assert isinstance(result.children[0], BandsPlotWidget)
