def test_outline_panel():
    """Test OutlinePanel class."""
    from aiidalab_qe.common.panel import OutlinePanel

    panel = OutlinePanel(identifier="test")
    assert panel.identifier == "test"
    parameters = panel.get_panel_value()
    assert parameters == {"test_run": False}
    parameters = {"test_run": True}
    panel.set_panel_value(parameters)
    assert panel.run.value is True
