import pytest


@pytest.mark.usefixtures("sssp")
def test_settings():
    """Test the settings of the xps app."""

    from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
    from ase.build import molecule
    from aiida.orm import StructureData

    configure_step = ConfigureQeAppWorkChainStep()
    # set the input structure
    h2o = molecule("H2O")
    h2o.center(vacuum=3.0)
    structure = StructureData(ase=h2o)
    configure_step.input_structure = structure
    # select xps
    configure_step.workchain_settings.properties["xps"].run.value = True
    # test get_panel_value
    configure_step.settings["xps"].structure_type.value = "molecule"
    # select the first elmement, which is O_1s
    configure_step.settings["xps"].core_level_list.children[0].value = True
    parameters = configure_step.settings["xps"].get_panel_value()
    print("parameters", parameters)
    assert parameters["structure_type"] == "molecule"
    assert parameters["core_level_list"] == ["O_1s"]
    assert (
        "not supported"
        in configure_step.settings["xps"].core_level_list.children[1].description
    )
    # set the parameters
    configure_step.settings["xps"].structure_type.value = "crystal"
    configure_step.settings["xps"].core_level_list.children[0].value = False
    configure_step.settings["xps"].set_panel_value(parameters)
    assert configure_step.settings["xps"].core_level_list.children[0].value is True
    assert configure_step.settings["xps"].structure_type.value == "molecule"
