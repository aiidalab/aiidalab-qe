import pytest

from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.configuration.model import ConfigurationModel
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION


def test_protocol():
    """Test the protocol.
    The protocol from workchain_settings will trigger the
    update of the protocol in advanced_settings.
    """
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    model.workchain.protocol = "fast"
    assert model.advanced.protocol == "fast"
    assert model.advanced.kpoints_distance == 0.5


def test_get_configuration_parameters():
    """Test the get_configuration_parameters method."""
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    parameters = config.get_configuration_parameters()
    parameters_ref = {
        "workchain": {
            **model.workchain.get_model_state(),
            "properties": model._get_properties(),
        },
        "advanced": model.advanced.get_model_state(),
    }
    assert parameters == parameters_ref


# TODO this test won't work unless the required pseudos are installed
# @pytest.mark.usefixtures("sssp")
# def test_set_configuration_parameters():
#     """Test the set_configuration_parameters method."""
#     model = ConfigurationModel()
#     config = ConfigureQeAppWorkChainStep(model=model)
#     config.render()
#     parameters = config.get_configuration_parameters()
#     parameters["workchain"]["relax_type"] = "positions"
#     parameters["advanced"]["pseudo_family"] = f"SSSP/{SSSP_VERSION}/PBE/efficiency"
#     config.set_configuration_parameters(parameters)
#     new_parameters = config.get_configuration_parameters()
#     assert parameters == new_parameters
#     # test pseudodojo
#     parameters["advanced"]["pseudo_family"] = (
#         f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf"
#     )
#     config.set_configuration_parameters(parameters)
#     new_parameters = config.get_configuration_parameters()
#     assert parameters == new_parameters


def test_panel():
    """Dynamic add/remove the panel based on the workchain settings."""
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    assert len(config.tab.children) == 2
    parameters = config.get_configuration_parameters()
    assert "bands" not in parameters
    model.get_model("bands").include = True
    assert len(config.tab.children) == 3
    parameters = config.get_configuration_parameters()
    assert "bands" in parameters


def test_reminder_info():
    """Dynamic add/remove the reminder text based on the workchain settings."""
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    outlines = config.workchain_settings.property_children[1:]
    bands_info = outlines[0].children[1]
    assert bands_info.value == ""
    bands_model = model.get_model("bands")
    bands_model.include = True
    assert bands_info.value == "Customize bands settings in the panel above if needed"
    bands_model.include = False
    assert bands_info.value == ""
