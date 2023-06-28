def test_result(generate_bands_workchain):
    from widget_bandsplot import BandsPlotWidget

    from aiidalab_qe.app.plugins.bands.result import Result, export_bands_data

    wkchain = generate_bands_workchain()
    #
    data = export_bands_data(wkchain.node)
    assert data is not None
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
    assert isinstance(result.children[0], BandsPlotWidget)
