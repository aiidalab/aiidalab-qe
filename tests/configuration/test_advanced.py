def test_smearing_settings():
    """Test SmearningSettings widget."""
    from aiidalab_qe.app.configuration.advanced import SmearingSettings

    # Check initial
    w = SmearingSettings()

    assert w.value.get("degauss") == 0.01
    assert w.value.get("smearing") == "cold"

    # Check values after changing protocol (currently the smearing not changed upon protocol)
    w.protocol = "fast"

    assert w.value.get("degauss") == 0.01
    assert w.value.get("smearing") == "cold"

    # Check changing value of sub-widgets changes the settings
    w.degauss.value = 0.03
    w.smearing.value = "methfessel-paxton"

    assert w.value.get("degauss") == 0.03
    assert w.value.get("smearing") == "methfessel-paxton"

    # Check reset
    w.reset()

    assert w.protocol == "moderate"
    assert w.value.get("degauss") == 0.01
    assert w.value.get("smearing") == "cold"


def test_advanced_kpoints_settings():
    """Test kpoint setting of advanced setting widget."""
    from aiidalab_qe.app.configuration.advanced import AdvancedSettings

    w = AdvancedSettings()

    # Check the disable of is bind to override switch
    assert w.kpoints_distance.disabled is True

    w.override.value = True
    assert w.kpoints_distance.disabled is False

    assert w.value.get("kpoints_distance") == 0.15

    w.protocol = "fast"

    assert w.value.get("kpoints_distance") == 0.5

    # Check changing value of sub-widgets changes the settings
    w.kpoints_distance.value = 0.22
    assert w.value.get("kpoints_distance") == 0.22

    # Check reset
    w.reset()

    assert w.protocol == "moderate"
    assert w.value.get("kpoints_distance") == 0.15


def test_tot_charge():
    """Test TotCharge widget."""
    from aiidalab_qe.app.configuration.advanced import TotalCharge

    w = TotalCharge()

    assert w.value == 0.0

    w.charge.value = 1.0
    assert w.value == 1.0

    # Check reset
    w.reset()

    assert w.value == 0.0
