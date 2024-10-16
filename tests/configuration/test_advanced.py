import pytest

from aiidalab_qe.app.main import App


def test_advanced_default():
    """Test default behavior of advanced setting."""
    app = App(qe_auto_setup=False)
    model = app.config_model

    # Test override functionality in advanced settings
    model.advanced.override = True
    model.workchain.protocol = "fast"
    model.advanced.smearing.type = "methfessel-paxton"
    model.advanced.smearing.degauss = 0.03
    model.advanced.kpoints_distance = 0.22

    # Reset values to default after removing override
    model.advanced.override = False

    assert model.advanced.smearing.type == "cold"
    assert model.advanced.smearing.degauss == 0.01
    assert model.advanced.kpoints_distance == 0.5


def test_advanced_smearing_settings():
    """Test Smearing Settings."""
    app = App(qe_auto_setup=False)
    model = app.config_model
    advanced = app.configure_step.advanced_settings
    advanced.smearing.render()

    # Test widget disable state on override
    assert advanced.smearing.degauss.disabled is True
    assert advanced.smearing.smearing.disabled is True

    model.advanced.override = True

    assert advanced.smearing.degauss.disabled is False
    assert advanced.smearing.smearing.disabled is False

    assert model.advanced.smearing.type == "cold"
    assert model.advanced.smearing.degauss == 0.01

    # Test protocol-dependent smearing change
    model.workchain.protocol = "fast"

    assert model.advanced.smearing.type == "cold"
    assert model.advanced.smearing.degauss == 0.01

    # Check reset
    model.advanced.smearing.type = "gaussian"
    model.advanced.smearing.degauss = 0.05
    model.advanced.smearing.reset()

    assert model.workchain.protocol == "fast"  # reset does not apply to protocol
    assert model.advanced.smearing.type == "cold"
    assert model.advanced.smearing.degauss == 0.01


def test_advanced_kpoints_settings():
    """Test kpoint setting of advanced setting widget."""
    app = App(qe_auto_setup=False)
    model = app.config_model
    advanced = app.configure_step.advanced_settings
    advanced.render()

    # Check the disable of is bind to override switch
    assert advanced.kpoints_distance.disabled is True

    model.advanced.override = True
    assert advanced.kpoints_distance.disabled is False

    assert model.advanced.kpoints_distance == 0.15

    model.workchain.protocol = "fast"
    assert model.advanced.kpoints_distance == 0.5

    # Check reset
    model.advanced.kpoints_distance = 0.1
    model.advanced.reset()

    assert model.workchain.protocol == "fast"  # reset does not apply to protocol
    assert model.advanced.kpoints_distance == 0.5


def test_advanced_tot_charge_settings():
    """Test TotCharge widget."""
    app = App(qe_auto_setup=False)
    model = app.config_model
    advanced = app.configure_step.advanced_settings
    advanced.render()

    # Check the disable of is bind to override switch
    assert advanced.total_charge.disabled is True

    model.advanced.override = True
    assert advanced.total_charge.disabled is False

    assert model.advanced.total_charge == 0.0

    # Check reset
    model.advanced.total_charge = 1.0
    model.advanced.reset()

    assert model.advanced.total_charge == 0.0


@pytest.mark.usefixtures("sssp")
def test_advanced_kpoints_mesh(generate_structure_data):
    """Test Mesh Grid HTML widget."""
    app = App(qe_auto_setup=False)
    model = app.config_model

    structure = generate_structure_data(name="silicon")
    model.input_structure = structure

    model.advanced.override = True
    assert model.advanced.mesh_grid == "Mesh [14, 14, 14]"

    # change protocol
    model.workchain.protocol = "fast"
    assert model.advanced.mesh_grid == "Mesh [5, 5, 5]"


@pytest.mark.usefixtures("sssp")
def test_advanced_hubbard_settings(generate_structure_data):
    """Test Hubbard widget."""
    app = App(qe_auto_setup=False)

    model = app.config_model
    hubbard_model = model.advanced.hubbard

    structure = generate_structure_data(name="LiCoO2")
    model.input_structure = structure

    # Activate Hubbard U widget
    hubbard_model.is_active = True
    assert hubbard_model.input_labels == ["Co - 3d", "O - 2p", "Li - 2s"]

    # Render the settings widget to allow for widget testing
    hubbard = app.configure_step.advanced_settings.hubbard
    hubbard.render()

    # Change the Hubbard U parameters for Co, O, and Li
    hubbard_parameters = hubbard.hubbard_widget.children[1:]  # type: ignore
    hubbard_parameters[0].value = 1  # Co - 3d
    hubbard_parameters[1].value = 2  # O - 2p
    hubbard_parameters[2].value = 3  # Li - 2s

    assert hubbard_model.parameters == {
        "Co - 3d": 1.0,
        "O - 2p": 2.0,
        "Li - 2s": 3.0,
    }

    # The widget hierarchy for eigenvalues:
    # - w.hubbard_widget.eigen_values_widget.children[0]: List of eigenvalues for Co
    # - w.hubbard_widget.eigen_values_widget.children[0].children[1]: Widgets for up and down spin
    # - w.hubbard_widget.eigen_values_widget.children[0].children[1].children[0]: Widget for up spin
    # - w.hubbard_widget.eigen_values_widget.children[0].children[1].children[0].children[1]: Widget for eigenvalue 1 (3d range: 1 to 5)

    # Check eigenvalues are empty
    # assert hubbard_model.eigenvalues == []  # TODO should they be?

    # Check there is only eigenvalues for Co (Transition metal)
    hubbard_model.eigenvalues_label = True

    assert len(hubbard_model.elements) == 1
    assert len(hubbard_model.eigenvalues) == 1

    Co_eigenvalues = hubbard.eigenvalues_widget.children[0].children[1]  # type: ignore
    Co_spin_down_row = Co_eigenvalues.children[1]
    Co_spin_down_row.children[1].value = "1"
    Co_spin_down_row.children[3].value = "1"
    Co_spin_down_row.children[5].value = "1"

    assert hubbard_model.get_active_eigenvalues() == [
        [1, 1, "Co", 1],
        [3, 1, "Co", 1],
        [5, 1, "Co", 1],
    ]
