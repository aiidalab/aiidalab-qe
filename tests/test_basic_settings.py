def test_get_panel_value():
    """Test get_panel_value."""
    from aiidalab_qe.app.configure.configure import BasicSettings

    wg = BasicSettings()
    parameters = wg.get_panel_value()
    parameters_ref = {
        "electronic_type": "metal",
        "spin_type": "none",
        "protocol": "moderate",
    }
    assert parameters == parameters_ref


def test_set_panel_value():
    """Test set_panel_value."""
    from aiidalab_qe.app.configure.configure import BasicSettings

    wg = BasicSettings()
    parameters_ref = {
        "electronic_type": "insulator",
        "spin_type": "collinear",
        "protocol": "fast",
    }
    wg.set_panel_value(parameters_ref)
    parameters = wg.get_panel_value()
    assert parameters == parameters_ref
