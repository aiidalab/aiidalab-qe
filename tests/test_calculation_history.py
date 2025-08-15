from aiidalab_qe.app.utils.calculation_history import CalculationHistory


def test_calculation_history(generate_qeapp_workchain):
    workchain = generate_qeapp_workchain()
    workchain.node.seal()

    calculation_history = CalculationHistory()
    calculation_history.load_table()
    assert len(calculation_history.table.data) >= 1
