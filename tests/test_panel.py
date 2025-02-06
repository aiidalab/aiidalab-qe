def test_panel_outline():
    """Test PanelOutline class."""
    from aiidalab_qe.common.panel import PluginOutline

    outline = PluginOutline()
    assert not outline.include.value
    outline.include.value = True
    assert outline.include.value
