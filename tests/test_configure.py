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
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    parameters = model.get_model_state()
    parameters_ref = {
        "workchain": {
            **model.workchain.get_model_state(),
            "properties": model._get_properties(),
        },
        "advanced": model.advanced.get_model_state(),
    }
    assert parameters == parameters_ref


@pytest.mark.usefixtures("aiida_profile_clean", "sssp", "pseudodojo")
def test_set_configuration_parameters():
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    parameters = model.get_model_state()
    parameters["workchain"]["relax_type"] = "positions"
    parameters["advanced"]["pseudo_family"] = f"SSSP/{SSSP_VERSION}/PBE/efficiency"
    model.set_model_state(parameters)
    new_parameters = model.get_model_state()
    assert parameters == new_parameters
    parameters["advanced"]["pseudo_family"] = (
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf"
    )
    model.set_model_state(parameters)
    new_parameters = model.get_model_state()
    assert parameters == new_parameters


def test_panel():
    """Dynamic add/remove the panel based on the workchain settings."""
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    assert len(config.tab.children) == 2
    parameters = model.get_model_state()
    assert "bands" not in parameters
    model.get_model("bands").include = True
    assert len(config.tab.children) == 3
    parameters = model.get_model_state()
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
