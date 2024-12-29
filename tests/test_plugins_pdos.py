def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bands_pdos import BandsPdosWidget
    from aiidalab_qe.plugins.electronic_structure.result import (
        ElectronicStructureResultsModel,
        ElectronicStructureResultsPanel,
    )

    workchain = generate_qeapp_workchain(run_bands=False)
    model = ElectronicStructureResultsModel()
    panel = ElectronicStructureResultsPanel(model=model)
    model.process_uuid = workchain.node.uuid

    assert model.title == "Electronic PDOS"
    assert model.identifiers == ["pdos"]

    panel.render()

    assert len(panel.results_container.children) == 1  # only pdos, so no controls

    widget = panel.bands_pdos_container.children[0]  # type: ignore
    model = widget._model

    assert isinstance(widget, BandsPdosWidget)
    assert isinstance(widget.plot, go.FigureWidget)

    assert not model.bands_data
    assert model.pdos_data

    assert widget.plot.layout.xaxis.title.text == "Density of states (eV)"
    assert widget.plot.layout.yaxis.title.text is None
