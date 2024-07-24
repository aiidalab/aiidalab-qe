import pytest


@pytest.mark.usefixtures("sssp")
def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go
    from aiidalab_qe.common.bandpdoswidget import BandPdosWidget
    from aiidalab_qe.plugins.bands.result import Result

    wkchain = generate_qeapp_workchain()
    # generate structure for scf calculation
    result = Result(wkchain.node)
    result._update_view()
    assert isinstance(result.children[0], BandPdosWidget)
    assert isinstance(result.children[0].bandsplot_widget, go.FigureWidget)

    # Check if data is correct
    assert result.children[0].bands_data is not None
    assert result.children[0].bands_data["pathlabels"] is not None
    assert result.children[0].pdos_data is None

    # Check Bands axis
    assert result.children[0].bandsplot_widget.layout.xaxis.title.text == "k-points"
    assert (
        result.children[0].bandsplot_widget.layout.yaxis.title.text
        == "Electronic Bands (eV)"
    )
    assert isinstance(
        result.children[0].bandsplot_widget.layout.xaxis.rangeslider,
        go.layout.xaxis.Rangeslider,
    )
    assert result.children[0].bands_data["pathlabels"][0] == list(
        result.children[0].bandsplot_widget.layout.xaxis.ticktext
    )


@pytest.mark.usefixtures("sssp")
def test_structure_1d(generate_qeapp_workchain, generate_structure_data):
    structure = generate_structure_data("silicon", pbc=(True, False, False))
    wkchain = generate_qeapp_workchain(structure=structure)
    assert "bands_kpoints_distance" not in wkchain.inputs.bands
    assert "bands_kpoints" in wkchain.inputs.bands
    assert len(wkchain.inputs.bands.bands_kpoints.labels) == 2


@pytest.mark.usefixtures("sssp")
def test_structure_2d(generate_qeapp_workchain, generate_structure_data):
    structure = generate_structure_data("silicon", pbc=(True, True, False))
    wkchain = generate_qeapp_workchain(structure=structure)
    assert "bands_kpoints_distance" not in wkchain.inputs.bands
    assert "bands_kpoints" in wkchain.inputs.bands
    assert len(wkchain.inputs.bands.bands_kpoints.labels) == 4
