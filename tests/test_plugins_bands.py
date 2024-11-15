def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandsPdosWidget
    from aiidalab_qe.plugins.bands.result import BandsResults, BandsResultsModel

    workchain = generate_qeapp_workchain()
    model = BandsResultsModel()
    model.process_node = workchain.node
    result = BandsResults(model=model)
    result.render()

    widget = result.children[0]
    model = widget._model

    assert isinstance(widget, BandsPdosWidget)
    assert isinstance(widget.plot, go.FigureWidget)

    # Check if data is correct
    assert not model.pdos_data
    assert model.bands_data
    assert model.bands_data["pathlabels"]  # type: ignore

    # Check Bands axis
    assert widget.plot.layout.xaxis.title.text == "k-points"
    assert widget.plot.layout.yaxis.title.text == "Electronic Bands (eV)"
    assert isinstance(widget.plot.layout.xaxis.rangeslider, go.layout.xaxis.Rangeslider)
    assert model.bands_data["pathlabels"][0] == list(widget.plot.layout.xaxis.ticktext)  # type: ignore


def test_structure_1d(generate_qeapp_workchain, generate_structure_data):
    structure = generate_structure_data("silicon", pbc=(True, False, False))
    workchain = generate_qeapp_workchain(structure=structure)
    assert "bands_kpoints_distance" not in workchain.inputs.bands.bands
    assert "bands_kpoints" in workchain.inputs.bands.bands
    assert len(workchain.inputs.bands.bands.bands_kpoints.labels) == 2
    assert workchain.inputs.bands.bands.bands_kpoints.labels == [(0, "Γ"), (9, "X")]


def test_structure_2d(generate_qeapp_workchain, generate_structure_data):
    structure = generate_structure_data("MoS2", pbc=(True, True, False))
    workchain = generate_qeapp_workchain(structure=structure)
    assert "bands_kpoints_distance" not in workchain.inputs.bands.bands
    assert "bands_kpoints" in workchain.inputs.bands.bands
    assert len(workchain.inputs.bands.bands.bands_kpoints.labels) == 4
    assert workchain.inputs.bands.bands.bands_kpoints.labels == [
        (0, "Γ"),
        (11, "M"),
        (18, "K"),
        (31, "Γ"),
    ]
