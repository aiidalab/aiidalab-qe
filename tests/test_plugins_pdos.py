import pytest


@pytest.mark.usefixtures("sssp")
def test_result(generate_qeapp_workchain):
    import plotly.graph_objects as go

    from aiidalab_qe.common.bandpdoswidget import BandPdosWidget
    from aiidalab_qe.plugins.pdos.result import Result

    wkchain = generate_qeapp_workchain()
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
    assert isinstance(result.children[0], BandPdosWidget)
    assert isinstance(result.children[0].bandsplot_widget, go.FigureWidget)

    # Check if data is correct
    assert result.children[0].bands_data is None
    assert result.children[0].pdos_data is not None

    # Check PDOS settings is not None

    # Check Bands axis
    assert (
        result.children[0].bandsplot_widget.layout.xaxis.title.text
        == "Density of states (eV)"
    )
    assert result.children[0].bandsplot_widget.layout.yaxis.title.text is None
