def test_protocol():
    """Test the protocol.
    The value should be the one in
    aiida_quantumespresso/workflows/protocols/pw/base.yaml
    """
    from aiidalab_qe.app.configure.advanced import AdvancedSettings

    wg = AdvancedSettings()
    assert wg.kpoints_distance.value == 0.15
    assert wg.degauss.value == 0.01
    assert wg.smearing.value == "cold"
    assert wg.tot_charge.value == 0
    # swich to fast protocol
    wg.workchain_protocol.value = "fast"
    assert wg.kpoints_distance.value == 0.5


def test_get_panel_value():
    """Test get_panel_value."""
    from aiidalab_qe.app.configure.advanced import AdvancedSettings

    wg = AdvancedSettings()
    parameters = wg.get_panel_value()
    parameters_ref = {
        "pseudo_family": "SSSP/1.2/PBEsol/efficiency",
        "kpoints_distance": 0.15,
        "pw": {
            "parameters": {
                "SYSTEM": {
                    "degauss": 0.01,
                    "smearing": "cold",
                    "tot_charge": 0,
                }
            },
        },
    }
    assert parameters == parameters_ref


def test_set_panel_value():
    """Test set_panel_value."""
    from aiidalab_qe.app.configure.advanced import AdvancedSettings

    wg = AdvancedSettings()
    parameters_ref = {
        "pseudo_family": "SSSP/1.2/PBE/precision",
        "kpoints_distance": 0.25,
        "pw": {
            "parameters": {
                "SYSTEM": {
                    "degauss": 0.02,
                    "smearing": "gaussian",
                    "tot_charge": 1,
                }
            },
        },
    }
    wg.set_panel_value(parameters_ref)
    parameters = wg.get_panel_value()
    assert parameters == parameters_ref
