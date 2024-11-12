def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandPdosWidget
    from aiidalab_qe.plugins.pdos.result import PdosResults, PdosResultsModel

    workchain = generate_qeapp_workchain()
    # generate structure for scf calculation
    model = PdosResultsModel()
    model.process_node = workchain.node
    result = PdosResults(model=model)
    result.render()

    widget = result.children[0]
    model = widget._model

    assert isinstance(widget, BandPdosWidget)

    widget.plot_button.click()
    assert isinstance(widget.plot, go.FigureWidget)

    # Check if data is correct
    assert not model.bands_data
    assert model.pdos_data

    # Check PDOS settings is not None

    # Check Bands axis
    assert widget.plot.layout.xaxis.title.text == "Density of states (eV)"
    assert widget.plot.layout.yaxis.title.text is None
