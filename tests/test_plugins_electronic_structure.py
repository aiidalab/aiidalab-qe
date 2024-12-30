def test_electronic_structure(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandsPdosWidget
    from aiidalab_qe.plugins.electronic_structure.result import (
        ElectronicStructureResultsModel,
        ElectronicStructureResultsPanel,
    )

    workchain = generate_qeapp_workchain()
    model = ElectronicStructureResultsModel()
    panel = ElectronicStructureResultsPanel(model=model)
    model.process_uuid = workchain.node.uuid

    assert model.title == "Electronic bands + PDOS"
    assert model.identifiers == ["bands", "pdos"]

    panel.render()

    assert len(panel.results_container.children) == 2  # has controls

    widget = panel.bands_pdos_container.children[0]  # type: ignore
    model = widget._model

    assert isinstance(widget, BandsPdosWidget)
    assert isinstance(widget.plot, go.FigureWidget)

    # Check if data is correct

    assert model.bands_data is not None
    assert model.bands_data["pathlabels"] is not None  # type: ignore
    assert model.pdos_data is not None

    # Check Bands axis
    assert widget.plot.layout.xaxis.title.text == "k-points"
    assert widget.plot.layout.xaxis2.title.text == "Density of states"
    assert widget.plot.layout.yaxis.title.text == "Electronic Bands (eV)"
    assert isinstance(
        widget.plot.layout.xaxis.rangeslider,
        go.layout.xaxis.Rangeslider,
    )
    assert model.bands_data["pathlabels"][0] == list(  # type: ignore
        widget.plot.layout.xaxis.ticktext
    )
