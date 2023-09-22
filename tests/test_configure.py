def test_protocol():
    """Test the protocol.
    The protocol from workchain_settings will trigger the
    update of the protocol in advanced_settings.
    """
    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    wg.workchain_settings.workchain_protocol.value = "fast"
    assert wg.advanced_settings.protocol == "fast"
    assert wg.advanced_settings.kpoints_distance.value == 0.5


def test_get_configuration_parameters():
    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    parameters = wg.get_configuration_parameters()
    parameters_ref = {
        "workchain": wg.workchain_settings.get_panel_value(),
        "advanced": wg.advanced_settings.get_panel_value(),
    }
    assert parameters == parameters_ref


def test_panel():
    """Dynamic add/remove the panel based on the workchain settings."""
    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    assert len(wg.tab.children) == 2
    parameters = wg.get_configuration_parameters()
    assert "bands" not in parameters
    wg.workchain_settings.properties["bands"].run.value = True
    assert len(wg.tab.children) == 3
    parameters = wg.get_configuration_parameters()
    assert "bands" in parameters
