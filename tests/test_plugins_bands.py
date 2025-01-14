def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandsPdosWidget
    from aiidalab_qe.plugins.electronic_structure.result import (
        ElectronicStructureResultsModel,
        ElectronicStructureResultsPanel,
    )

    workchain = generate_qeapp_workchain(run_pdos=False)
    model = ElectronicStructureResultsModel()
    panel = ElectronicStructureResultsPanel(model=model)

    model.process_uuid = workchain.node.uuid

    assert model.title == "Electronic bands"
    assert model.identifiers == ["bands"]

    panel.render()

    assert len(panel.results_container.children) == 1  # only bands, so no controls

    widget = panel.bands_pdos_container.children[0]  # type: ignore
    model = widget._model

    assert isinstance(widget, BandsPdosWidget)
    assert isinstance(widget.plot, go.FigureWidget)

    assert not model.pdos_data
    assert model.bands_data
    assert model.bands_data["pathlabels"]

    assert widget.plot.layout.xaxis.title.text == "k-points"
    assert widget.plot.layout.yaxis.title.text == "Electronic Bands (eV)"
    assert isinstance(widget.plot.layout.xaxis.rangeslider, go.layout.xaxis.Rangeslider)
    assert model.bands_data["pathlabels"][0] == list(widget.plot.layout.xaxis.ticktext)


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
