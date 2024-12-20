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
    assert smearing_model.degauss == 0.01
    assert model.kpoints_distance == 0.5


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
    assert smearing_model.degauss == 0.01

    # Check reset
    smearing_model.type = "gaussian"
    smearing_model.degauss = 0.05

    smearing_model.reset()

    assert smearing_model.type == "cold"
    assert smearing_model.degauss == 0.01


def test_advanced_kpoints_settings():
    """Test kpoints setting of advanced setting widget."""
    model = AdvancedConfigurationSettingsModel()
    _ = AdvancedConfigurationSettingsPanel(model=model)

    model.protocol = "moderate"
    assert model.kpoints_distance == 0.15

    model.protocol = "fast"
    assert model.kpoints_distance == 0.5

    # Check reset
    model.kpoints_distance = 0.1
    model.reset()

    assert model.protocol == "fast"  # reset does not apply to protocol
    assert model.kpoints_distance == 0.5


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
    assert model.mesh_grid == "Mesh [5, 5, 5]"


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
    hubbard_parameters = hubbard.hubbard_widget.children[1:]  # type: ignore
    hubbard_parameters[0].value = 1  # Co - 3d
    hubbard_parameters[1].value = 2  # O - 2p
    hubbard_parameters[2].value = 3  # Li - 2s

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

    assert model.get_active_eigenvalues() == [
        [1, 1, "Co", 1],
        [3, 1, "Co", 1],
        [5, 1, "Co", 1],
    ]
