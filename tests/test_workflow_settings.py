def test_get_panel_value():
    """Test get_panel_value."""
    from aiidalab_qe.app.configure.workflow import WorkChainSettings

    wg = WorkChainSettings()
    parameters = wg.get_panel_value()
    parameters_ref = {"relax_type": "positions_cell", "properties": {}}
    assert parameters == parameters_ref


def test_set_panel_value():
    """Test set_panel_value."""
    from aiidalab_qe.app.configure.workflow import WorkChainSettings

    wg = WorkChainSettings()
    parameters_ref = {"relax_type": "positions", "properties": {}}
    wg.set_panel_value(parameters_ref)
    parameters = wg.get_panel_value()
    assert parameters == parameters_ref
