def test_result(generate_qeapp_workchain):
    from widget_bandsplot import BandsPlotWidget

    from aiidalab_qe.plugins.bands.result import Result, export_bands_data

    wkchain = generate_qeapp_workchain()
    data = export_bands_data(wkchain.node.outputs.bands)
    assert data is not None
    # generate structure for scf calculation
    result = Result(wkchain.node)
    assert result.identifier == "bands"
    result._update_view()
    assert isinstance(result.children[0], BandsPlotWidget)


def test_structure_1d(generate_qeapp_workchain, generate_structure_data):
    structure = generate_structure_data("silicon", pbc=(True, False, False))
    wkchain = generate_qeapp_workchain(structure=structure)
    assert "bands_kpoints_distance" not in wkchain.inputs.bands
    assert "bands_kpoints" in wkchain.inputs.bands
    assert len(wkchain.inputs.bands.bands_kpoints.labels) == 2


def test_structure_2d(generate_qeapp_workchain, generate_structure_data):
    structure = generate_structure_data("silicon", pbc=(True, True, False))
    wkchain = generate_qeapp_workchain(structure=structure)
    assert "bands_kpoints_distance" not in wkchain.inputs.bands
    assert "bands_kpoints" in wkchain.inputs.bands
    assert len(wkchain.inputs.bands.bands_kpoints.labels) == 4
