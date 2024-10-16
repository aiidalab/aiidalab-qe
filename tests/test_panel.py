def test_panel_outline():
    """Test PanelOutline class."""
    from aiidalab_qe.common.panel import PanelOutline

    outline = PanelOutline(identifier="test")
    assert outline.identifier == "test"
    assert not outline.include.value
    outline.include.value = True
    assert outline.include.value
