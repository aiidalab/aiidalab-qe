def test_panel_class():
    """Test Panel class."""
    from aiidalab_qe.app.common.panel import Panel

    panel = Panel(identifier="test")
    assert panel.identifier == "test"
    parameters = {"test": "test"}
    panel.set_panel_value(parameters)
    parameters = panel.get_panel_value()


def test_outline_panel():
    """Test OutlinePanel class."""
    from aiidalab_qe.app.common.panel import OutlinePanel

    panel = OutlinePanel(identifier="test")
    assert panel.identifier == "test"
    parameters = panel.get_panel_value()
    assert parameters == {"test_run": False}
    parameters = {"test_run": True}
    panel.set_panel_value(parameters)
    assert panel.run.value is True


def test_result_panel():
    """Test ResultPanel class.""" ""
    from aiidalab_qe.app.common.panel import ResultPanel

    panel = ResultPanel(identifier="test")
    assert panel.identifier == "test"
    panel._update_view()
