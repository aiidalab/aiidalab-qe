def test_result(generate_pdos_workchain):
    from aiidalab_qe.app.plugins.pdos.result import Result

    wkchain = generate_pdos_workchain()
    # generate structure for scf calculation
    result = Result(node=wkchain.node)
    result._update_view()
