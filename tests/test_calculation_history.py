def test_calculation_history(sssp, generate_qeapp_workchain):
    from aiidalab_qe.app.utils.calculation_history import CalculationHistory

    workchain = generate_qeapp_workchain()
    workchain.node.seal()

    calculation_history = CalculationHistory()
    calculation_history.load_table()
    assert len(calculation_history.table.data) >= 1
