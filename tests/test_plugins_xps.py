from aiidalab_qe.app.configuration.model import ConfigurationStepModel


def test_settings():
    """Test the settings of the xps app."""

    from ase.build import molecule

    from aiida.orm import StructureData
    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep

    # TODO use submit_app_generator(properties=["xps"]) instead? See xas test
    config_model = ConfigurationStepModel()
    _ = ConfigureQeAppWorkChainStep(model=config_model)

    # Set the input structure
    h2o = molecule("H2O")
    h2o.center(vacuum=3.0)
    structure = StructureData(ase=h2o)
    config_model.input_structure = structure

    # Select xps
    xps_model = config_model.get_model("xps")
    xps_model.include = True
    xps_model.update()

    # Test getting the model state
    xps_model.structure_type = "molecule"
    xps_model.core_levels["O_1s"] = True
    parameters = xps_model.get_model_state()
    assert parameters["structure_type"] == "molecule"
    assert parameters["core_level_list"] == ["O_1s"]
    assert "H" not in xps_model.get_supported_core_levels()  # H not supported

    # Test setting the model state
    xps_model.structure_type = "crystal"
    xps_model.core_levels["O_1s"] = False
    xps_model.set_model_state(parameters)
    assert xps_model.core_levels["O_1s"] is True
    assert xps_model.structure_type == "molecule"
