def test_result(generate_pdos_workchain):
    from aiidalab_qe.app.plugins.pdos.result import Result, export_pdos_data

    wkchain = generate_pdos_workchain()
    data = export_pdos_data(wkchain.node)
    assert data is not None
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
    assert len(result.children) == 2


def test_result_spin(generate_pdos_workchain):
    from aiidalab_qe.app.plugins.pdos.result import Result, export_pdos_data

    wkchain = generate_pdos_workchain(spin=True)
    data = export_pdos_data(wkchain.node)
    assert data is not None
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
    assert len(result.children) == 2


def test_result_group_by(generate_pdos_workchain):
    from aiidalab_qe.app.plugins.pdos.result import Result, export_pdos_data

    wkchain = generate_pdos_workchain()
    data = export_pdos_data(wkchain.node)
    assert data is not None
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
    result.children[0].children[0].children[1].value = "angular"
