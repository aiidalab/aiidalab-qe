from aiida.common.extendeddicts import AttributeDict


def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget

    workchain = generate_qeapp_workchain()

    # NOTE the actual widget fails because the workchain is not actually attached
    # to the QeAppWorkchain, so the bands widget receives `None` and raises an
    # exception. Instead, we mock the render behavior, but bypass the node fetching
    # by setting the node directly from the outputs of the generated workchain.
    # TODO rethink test

    pdos_node = AttributeDict(workchain.outputs["pdos"])
    model = BandsPdosModel()
    widget = BandsPdosWidget(model=model, pdos=pdos_node)
    widget.render()

    assert isinstance(widget, BandsPdosWidget)
    assert isinstance(widget.plot, go.FigureWidget)

    # Check if data is correct
    assert not model.bands_data
    assert model.pdos_data

    # Check PDOS settings is not None

    # Check Bands axis
    assert widget.plot.layout.xaxis.title.text == "Density of states (eV)"
    assert widget.plot.layout.yaxis.title.text is None
