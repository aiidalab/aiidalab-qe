from aiida.common.extendeddicts import AttributeDict


def test_electronic_structure(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget

    workchain = generate_qeapp_workchain()

    # NOTE the actual widget fails because the workchain is not actually attached
    # to the QeAppWorkchain, so the bands widget receives `None` and raises an
    # exception. Instead, we mock the render behavior, but bypass the node fetching
    # by setting the node directly from the outputs of the generated workchain.
    # TODO rethink test

    bands_node = workchain.outputs["bands"]["bands"]
    pdos_node = AttributeDict(workchain.outputs["pdos"])
    model = BandsPdosModel()
    widget = BandsPdosWidget(model=model, bands=bands_node, pdos=pdos_node)
    widget.render()

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
