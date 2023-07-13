def test_protocol():
    """Test the protocol.
    The protocol from basic_settings will trigger the
    update of the protocol in advance_settings.
    """
    from aiidalab_qe.app.configure.configure import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    wg.basic_settings.workchain_protocol.value = "fast"
    assert wg.advance_settings.workchain_protocol.value == "fast"
    assert wg.advance_settings.kpoints_distance.value == 0.5


def test_get_input_parameters():
    from aiidalab_qe.app.configure.configure import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    parameters = wg.get_input_parameters()
    parameters_ref = {
        "basic": wg.basic_settings.get_panel_value(),
        "advanced": wg.advance_settings.get_panel_value(),
        "workflow": wg.workchain_settings.get_panel_value(),
    }
    assert parameters == parameters_ref


def test_set_input_parameters():
    """insulator, non-magnetic, no smearing
    the occupation type is set to fixed, smearing and degauss should not be set"""
    from aiidalab_qe.app.configure.configure import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    parameters_ref = {
        "basic": {
            "electronic_type": "insulator",
            "spin_type": "collinear",
            "protocol": "fast",
        },
        "advanced": {
            "pseudo_family": "SSSP/1.2/PBE/precision",
            "kpoints_distance": 0.25,
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "degauss": 0.02,
                        "smearing": "gaussian",
                        "tot_charge": 0,
                    }
                },
            },
        },
        "workflow": {
            "relax_type": "positions",
            "properties": {
                "bands": False,
                "pdos": False,
            },
        },
    }
    wg.set_input_parameters(parameters_ref)
    parameters = wg.get_input_parameters()
    assert parameters == parameters_ref


def test_panel():
    """Dynamic add/remove the panel based on the workchain settings."""
    from aiidalab_qe.app.configure.configure import ConfigureQeAppWorkChainStep

    wg = ConfigureQeAppWorkChainStep()
    assert len(wg.tab.children) == 3
    parameters = wg.get_input_parameters()
    assert "bands" not in parameters
    wg.workchain_settings.properties["bands"].run.value = True
    assert len(wg.tab.children) == 4
    parameters = wg.get_input_parameters()
    assert "bands" in parameters
