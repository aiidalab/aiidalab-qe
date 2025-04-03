from aiidalab_qe.app.configuration.advanced import (
    AdvancedConfigurationSettingsModel,
    AdvancedConfigurationSettingsPanel,
)


def test_advanced_default():
    """Test default behavior of advanced setting."""
    model = AdvancedConfigurationSettingsModel()
    _ = AdvancedConfigurationSettingsPanel(model=model)

    smearing_model = model.get_model("smearing")

    model.protocol = "fast"
    smearing_model.type = "methfessel-paxton"
    smearing_model.degauss = 0.03
    model.kpoints_distance = 0.22

    # Reset values to default w.r.t protocol
    model.reset()
    smearing_model.reset()

    assert smearing_model.type == "cold"
    assert smearing_model.degauss == 0.0275
    assert model.kpoints_distance == 0.3


def test_advanced_smearing_settings():
    """Test Smearing Settings."""

    model = AdvancedConfigurationSettingsModel()
    _ = AdvancedConfigurationSettingsPanel(model=model)

    smearing_model = model.get_model("smearing")

    assert smearing_model.type == "cold"
    assert smearing_model.degauss == 0.01

    # Test protocol-dependent smearing change
    model.protocol = "fast"

    assert smearing_model.type == "cold"
    assert smearing_model.degauss == 0.0275

    # Check reset
    smearing_model.type = "gaussian"
    smearing_model.degauss = 0.05

    smearing_model.reset()

    assert smearing_model.type == "cold"
    assert smearing_model.degauss == 0.0275


def test_advanced_kpoints_settings():
    """Test kpoints setting of advanced setting widget."""
    model = AdvancedConfigurationSettingsModel()
    _ = AdvancedConfigurationSettingsPanel(model=model)

    model.protocol = "balanced"
    assert model.kpoints_distance == 0.15

    model.protocol = "fast"
    assert model.kpoints_distance == 0.3

    # Check reset
    model.kpoints_distance = 0.1
    model.reset()

    assert model.protocol == "fast"  # reset does not apply to protocol
    assert model.kpoints_distance == 0.3


def test_advanced_molecule_settings(generate_structure_data):
    """Test kpoints setting of advanced setting widget."""
    model = AdvancedConfigurationSettingsModel()
    _ = AdvancedConfigurationSettingsPanel(model=model)

    # Create molecule
    structure = generate_structure_data(name="H2O", pbc=(False, False, False))
    model.input_structure = structure

    # Confirm the value of kpoints_distance is fixed
    assert model.kpoints_distance == 100.0

    model.protocol = "fast"
    assert model.kpoints_distance == 100.0

    # Check that reset is done w.r.t the molecule structure
    model.reset()
    assert model.protocol == "fast"  # reset does not apply to protocol
    assert model.kpoints_distance == 100


def test_advanced_tot_charge_settings():
    """Test TotCharge widget."""
    model = AdvancedConfigurationSettingsModel()
    _ = AdvancedConfigurationSettingsPanel(model=model)

    assert model.total_charge == 0.0

    # Check reset
    model.total_charge = 1.0
    model.reset()

    assert model.total_charge == 0.0


def test_advanced_kpoints_mesh(generate_structure_data):
    """Test Mesh Grid HTML widget."""
    model = AdvancedConfigurationSettingsModel()
    _ = AdvancedConfigurationSettingsPanel(model=model)

    structure = generate_structure_data(name="silicon")
    model.input_structure = structure

    assert model.mesh_grid == "Mesh [14, 14, 14]"

    # change protocol
    model.protocol = "fast"
    assert model.mesh_grid == "Mesh [7, 7, 7]"


def test_advanced_hubbard_settings(generate_structure_data):
    """Test Hubbard widget."""
    from aiidalab_qe.app.configuration.advanced.hubbard import (
        HubbardConfigurationSettingsModel,
        HubbardConfigurationSettingsPanel,
    )

    model = HubbardConfigurationSettingsModel()
    hubbard = HubbardConfigurationSettingsPanel(model=model)
    hubbard.render()

    structure = generate_structure_data(name="LiCoO2")
    model.input_structure = structure

    # Activate Hubbard U widget
    model.is_active = True
    assert model.orbital_labels == ["Co - 3d", "O - 2p", "Li - 2s"]

    # Change the Hubbard U parameters for Co, O, and Li
    hubbard_parameters = hubbard.hubbard_widget.children[:]  # type: ignore
    hubbard_parameters[0].children[0].value = 1  # Co - 3d
    hubbard_parameters[1].children[0].value = 2  # O - 2p
    hubbard_parameters[2].children[0].value = 3  # Li - 2s

    assert model.parameters == {
        "Co - 3d": 1.0,
        "O - 2p": 2.0,
        "Li - 2s": 3.0,
    }

    # The widget hierarchy for eigenvalues:
    # - hubbard.eigenvalues_widget.children[0]: List of eigenvalues for Co
    # - hubbard.eigenvalues_widget.children[0].children[1]: Widgets for up and down spin
    # - hubbard.eigenvalues_widget.children[0].children[1].children[0]: Widget for up spin
    # - hubbard.eigenvalues_widget.children[0].children[1].children[0].children[1]: Widget for eigenvalue 1 (3d range: 1 to 5)

    # Check eigenvalues are empty
    # assert model.hubbard.eigenvalues == []  # TODO should they be?

    # Check there is only eigenvalues for Co (Transition metal)
    model.has_eigenvalues = True
    assert len(model.applicable_kind_names) == 1
    assert len(model.eigenvalues) == 1

    Co_eigenvalues = hubbard.eigenvalues_widget.children[0].children[1]  # type: ignore
    Co_spin_down_row = Co_eigenvalues.children[1]
    Co_spin_down_row.children[1].value = "1"
    Co_spin_down_row.children[3].value = "1"
    Co_spin_down_row.children[5].value = "1"

    Co_spin_up_row = Co_eigenvalues.children[0]
    Co_spin_up_row.children[1].value = "0"
    Co_spin_up_row.children[3].value = "0"
    Co_spin_up_row.children[5].value = "0"
    assert model.get_active_eigenvalues() == [
        (1, 1, "Co", 0.0),
        (3, 1, "Co", 0.0),
        (5, 1, "Co", 0.0),
        (1, 2, "Co", 1.0),
        (3, 2, "Co", 1.0),
        (5, 2, "Co", 1.0),
    ]


def test_advanced_magnetic_settings(generate_structure_data):
    """Test Magnetization widget."""
    from aiida.orm import StructureData
    from aiidalab_qe.app.configuration.advanced.magnetization import (
        MagnetizationConfigurationSettingsModel,
        MagnetizationConfigurationSettingsPanel,
    )

    model = MagnetizationConfigurationSettingsModel()
    magnetic = MagnetizationConfigurationSettingsPanel(model=model)
    magnetic.render()

    model.family = "SSSP/1.3/PBE/efficiency"

    structure = generate_structure_data(name="LiCoO2")
    model.input_structure = structure
    model._update_default_moments()

    # The sssp fixture sets the number of valence electrons to 4 for all electrons. Here, we are only
    # testing that the correct logic is applied, i.e., 0.1 * number of electrons for elements with
    # default magnetic moment 0.
    assert model._defaults["moments"] == {"Li": 0.4, "Co": 5, "O": 0.4}

    structure = StructureData(
        cell=[
            [3.84737, 0.0, 0.0],
            [1.923685, 3.331920, 0.0],
            [1.923685, 1.110640, 3.141364],
        ]
    )
    structure.append_atom(position=(0.0, 0.0, 0.0), symbols="Ni", name="Ni1")
    structure.append_atom(
        position=(1.923685, 1.110640, 0.785341), symbols="Ni", name="Ni2"
    )
    structure.append_atom(position=(1.923685, 0.0, 2.356204), symbols="O", name="O1")
    structure.append_atom(position=(1.923685, 0.0, 0.785341), symbols="O", name="O2")

    model.input_structure = structure
    model._update_default_moments()

    assert model._defaults["moments"] == {"O1": 0.4, "O2": 0.4, "Ni1": 5, "Ni2": 5}
