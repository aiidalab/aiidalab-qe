def test_advanced_default():
    """Test default behavior of advanced setting widget."""
    from aiidalab_qe.app.configuration.advanced import AdvancedSettings

    w = AdvancedSettings()

    # Test override values and set back to default after un-tick the override checkbox
    w.override.value = True
    w.protocol = "fast"
    w.smearing.smearing.value = "methfessel-paxton"
    w.smearing.degauss.value = 0.03
    w.kpoints_distance.value = 0.22

    assert w.value.get("degauss") == 0.03
    assert w.value.get("smearing") == "methfessel-paxton"
    assert w.value.get("kpoints_distance") == 0.22

    w.override.value = False

    assert w.value.get("degauss") == 0.01
    assert w.value.get("smearing") == "cold"
    assert w.value.get("kpoints_distance") == 0.5


def test_advanced_smearing_settings():
    """Test SmearningSettings widget."""
    from aiidalab_qe.app.configuration.advanced import AdvancedSettings

    w = AdvancedSettings()

    # Check the disable of is bind to override switch of both degauss and smearing
    assert w.smearing.disabled is True
    assert w.smearing.degauss.disabled is True
    assert w.smearing.smearing.disabled is True

    w.override.value = True

    assert w.smearing.disabled is False
    assert w.smearing.degauss.disabled is False
    assert w.smearing.smearing.disabled is False

    assert w.value.get("degauss") == 0.01
    assert w.value.get("smearing") == "cold"

    # Check values after changing protocol (currently the smearing not changed upon protocol)
    w.protocol = "fast"

    assert w.value.get("degauss") == 0.01
    assert w.value.get("smearing") == "cold"

    # Check changing value of sub-widgets changes the settings
    # By set the widget value directly
    w.smearing.degauss.value = 0.03
    w.smearing.smearing.value = "methfessel-paxton"

    assert w.value.get("degauss") == 0.03
    assert w.value.get("smearing") == "methfessel-paxton"

    # Check reset
    w.reset()

    # the protocol will not be reset
    assert w.protocol == "fast"
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

    # the protocol will not be reset
    assert w.protocol == "fast"
    assert w.value.get("kpoints_distance") == 0.5


def test_advanced_tot_charge_settings():
    """Test TotCharge widget."""
    from aiidalab_qe.app.configuration.advanced import AdvancedSettings

    w = AdvancedSettings()

    # Check the disable of is bind to override switch
    assert w.total_charge.disabled is True

    w.override.value = True
    assert w.total_charge.disabled is False

    assert w.value.get("total_charge") == 0.0

    w.total_charge.value = 1.0
    assert w.value.get("total_charge") == 1.0

    # Check reset
    w.reset()

    assert w.value.get("total_charge") == 0.0


def test_advanced_kpoints_mesh():
    """Test Mesh Grid HTML widget."""
    from aiida import orm

    from aiidalab_qe.app.configuration.advanced import AdvancedSettings

    w = AdvancedSettings()

    # Create a StructureData for AdvancedSettings (Silicon)

    structure = orm.StructureData(
        cell=[
            [3.8401979337, 0.0000000000, 0.0000000000],
            [1.9200989668, 3.3257101909, 0.0000000000],
            [0.0000000000, -2.2171384943, 3.1355090603],
        ],
        pbc=[True, True, True],
    )
    structure.append_atom(
        position=[0.0000000000, 0.0000000000, 0.0000000000], symbols="Si"
    )
    structure.append_atom(
        position=[2.5601312956, 1.6544735986, 1.5677545302], symbols="Si"
    )

    w.input_structure = structure

    w.override.value = True
    assert w.mesh_grid.value == "Mesh [16, 16, 16]"

    # change protocol
    w.protocol = "fast"
<<<<<<< HEAD
    assert w.mesh_grid.value == "Mesh [5, 5, 5]"


def test_advanced_hubbard_widget():
    """Test Hubbard widget."""
    from aiida import orm

    from aiidalab_qe.app.configuration.advanced import AdvancedSettings

    w = AdvancedSettings()

    # Create a StructureData for AdvancedSettings (LiCoO2)

    a, b, c, d = 1.4060463552647, 0.81178124180108, 4.6012019181836, 1.6235624832021
    cell = [[a, -b, c], [0.0, d, c], [-a, -b, c]]
    sites = [
        ["Co", "Co", (0, 0, 0)],
        ["O", "O", (0, 0, 3.6020728736387)],
        ["O", "O", (0, 0, 10.201532881212)],
        ["Li", "Li", (0, 0, 6.9018028772754)],
    ]
    structure = orm.StructureData(cell=cell)

    for site in sites:
        structure.append_atom(position=site[2], symbols=site[0], name=site[1])

    w.input_structure = structure

    # Activate Hubbard U widget
    w.hubbard_widget.hubbard.value = True

    assert w.hubbard_widget.input_labels == ["Co - 3d", "O - 2p", "Li - 2s"]

    # Change the value of the Hubbard U for Co, O and Li
    w.hubbard_widget.hubbard_widget.children[1].children[0].value = 1
    w.hubbard_widget.hubbard_widget.children[2].children[0].value = 2
    w.hubbard_widget.hubbard_widget.children[3].children[0].value = 3

    assert w.hubbard_widget.hubbard_dict == {
        "hubbard_u": {"Co - 3d": 1.0, "O - 2p": 2.0, "Li - 2s": 3.0}
    }

    # Check eigenvalues are empty
    assert w.hubbard_widget.eigenvalues_dict == {}

    w.hubbard_widget.eigenvalues_label.value = True

    # Check there is only eigenvalues for Co (Transition metal)

    assert len(w.hubbard_widget.eigen_values_widget.children) == 1

    # The widget hierarchy for eigenvalues:
    # - w.hubbard_widget.eigen_values_widget.children[0]: List of eigenvalues for Co
    # - w.hubbard_widget.eigen_values_widget.children[0].children[1]: Widgets for up and down spin
    # - w.hubbard_widget.eigen_values_widget.children[0].children[1].children[0]: Widget for up spin
    # - w.hubbard_widget.eigen_values_widget.children[0].children[1].children[0].children[1]: Widget for eigenvalue 1 (3d range: 1 to 5)

    w.hubbard_widget.eigen_values_widget.children[0].children[1].children[0].children[
        1
    ].value = "1"
    w.hubbard_widget.eigen_values_widget.children[0].children[1].children[0].children[
        3
    ].value = "1"
    w.hubbard_widget.eigen_values_widget.children[0].children[1].children[0].children[
        5
    ].value = "1"

    assert w.hubbard_widget.eigenvalues_dict == {
        "starting_ns_eigenvalue": [[1, 1, "Co", 1], [3, 1, "Co", 1], [5, 1, "Co", 1]]
    }
=======
    assert w.mesh_grid.value == "Mesh [5, 5, 5]"
>>>>>>> a17893e786a289341c94ad69e53e6905201d5899
