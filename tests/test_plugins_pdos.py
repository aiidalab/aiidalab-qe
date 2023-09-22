def test_result(generate_qeapp_workchain):
    from aiidalab_qe.plugins.pdos.result import Result, export_pdos_data

    wkchain = generate_qeapp_workchain()
    data = export_pdos_data(wkchain.node.outputs.pdos)
    assert data is not None
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    assert result.identifier == "pdos"
    result._update_view()
    assert len(result.children) == 2


def test_result_spin(generate_qeapp_workchain):
    from aiidalab_qe.plugins.pdos.result import Result, export_pdos_data

    wkchain = generate_qeapp_workchain(spin_type="collinear")
    data = export_pdos_data(wkchain.node.outputs.pdos)
    assert data is not None
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
    assert len(result.children) == 2


def test_result_group_by(generate_qeapp_workchain):
    from aiidalab_qe.plugins.pdos.result import Result, export_pdos_data

    wkchain = generate_qeapp_workchain()
    data = export_pdos_data(wkchain.node.outputs.pdos)
    assert data is not None
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
    result.children[0].children[0].children[1].value = "angular"
