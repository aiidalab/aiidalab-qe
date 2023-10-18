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
    """Test the get_configuration_parameters method."""
    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    parameters = wg.get_configuration_parameters()
    parameters_ref = {
        "workchain": wg.workchain_settings.get_panel_value(),
        "advanced": wg.advanced_settings.get_panel_value(),
    }
    assert parameters == parameters_ref


def test_set_configuration_parameters():
    """Test the set_configuration_parameters method."""
    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    parameters = wg.get_configuration_parameters()
    parameters["workchain"]["relax_type"] = "positions"
    parameters["advanced"]["pseudo_family"] = "SSSP/1.2/PBE/efficiency"
    wg.set_configuration_parameters(parameters)
    new_parameters = wg.get_configuration_parameters()
    assert parameters == new_parameters
    # test pseudodojo
    parameters["advanced"]["pseudo_family"] = "PseudoDojo/0.4/PBEsol/SR/standard/upf"
    wg.set_configuration_parameters(parameters)
    new_parameters = wg.get_configuration_parameters()
    assert parameters == new_parameters


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
