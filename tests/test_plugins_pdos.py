def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandsPdosWidget
    from aiidalab_qe.plugins.pdos.result import PdosResultsModel, PdosResultsPanel

    workchain = generate_qeapp_workchain()
    model = PdosResultsModel()
    model.process_uuid = workchain.node.uuid
    result = PdosResultsPanel(model=model)
    result.render()

    widget = result.children[0]
    model = widget._model

    assert isinstance(widget, BandsPdosWidget)
    assert isinstance(widget.plot, go.FigureWidget)

    # Check if data is correct
    assert not model.bands_data
    assert model.pdos_data

    # Check PDOS settings is not None

    # Check Bands axis
    assert widget.plot.layout.xaxis.title.text == "Density of states (eV)"
    assert widget.plot.layout.yaxis.title.text is None
