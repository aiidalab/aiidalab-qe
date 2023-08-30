def test_smearing_settings():
    """Test SmearningSettings widget."""
    from aiidalab_qe.app.configuration.advanced import SmearingSettings

    # Check initial
    w = SmearingSettings()

    assert w.settings.get("degauss") == 0.01
    assert w.settings.get("smearing") == "cold"

    # Check values after changing protocol (currently the smearing not changed upon protocol)
    w.protocol = "fast"

    assert w.settings.get("degauss") == 0.01
    assert w.settings.get("smearing") == "cold"

    # Check changing value of sub-widgets changes the settings
    w.degauss.value = 0.03
    w.smearing.value = "methfessel-paxton"

    assert w.settings.get("degauss") == 0.03
    assert w.settings.get("smearing") == "methfessel-paxton"

    # Check reset
    w.reset()

    assert w.protocol == "moderate"
    assert w.settings.get("degauss") == 0.01
    assert w.settings.get("smearing") == "cold"


def test_kpoints_settings():
    """Test KpointSettings widget."""
    from aiidalab_qe.app.configuration.advanced import KpointSettings

    w = KpointSettings()

    assert w.settings.get("kpoints_distance") == 0.15

    w.protocol = "fast"

    assert w.settings.get("kpoints_distance") == 0.5

    # Check changing value of sub-widgets changes the settings
    w.distance.value = 0.22
    assert w.settings.get("kpoints_distance") == 0.22

    # Check reset
    w.reset()

    assert w.protocol == "moderate"
    assert w.settings.get("kpoints_distance") == 0.15
