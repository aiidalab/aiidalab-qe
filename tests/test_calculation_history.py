import pytest

@pytest.mark.usefixtures("aiida_profile_clean")
def test_calculation_history(generate_qeapp_workchain):
    from aiidalab_qe.app.utils.search_jobs import CalculationHistory

    workchain = generate_qeapp_workchain()
    workchain.node.seal()

    calculation_history = CalculationHistory()
    calculation_history.load_table()
    assert len(calculation_history.table.data) == 1