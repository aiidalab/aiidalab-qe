from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.configuration.model import ConfigurationModel
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION


def test_protocol():
    model = ConfigurationModel()
    _ = ConfigureQeAppWorkChainStep(model=model)
    workchain_model = model.get_model("workchain")
    advanced_model = model.get_model("advanced")
    workchain_model.protocol = "fast"
    assert advanced_model.protocol == "fast"
    assert advanced_model.kpoints_distance == 0.5


def test_get_configuration_parameters():
    model = ConfigurationModel()
    _ = ConfigureQeAppWorkChainStep(model=model)
    parameters = model.get_model_state()
    parameters_ref = {
        "workchain": {
            **model.get_model("workchain").get_model_state(),
            "relax_type": model.relax_type,
            "properties": model._get_properties(),
        },
        "advanced": model.get_model("advanced").get_model_state(),
    }
    assert parameters == parameters_ref


def test_set_configuration_parameters():
    model = ConfigurationModel()
    _ = ConfigureQeAppWorkChainStep(model=model)
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
    assert len(config.tabs.children) == 2
    parameters = model.get_model_state()
    assert "bands" not in parameters
    model.get_model("bands").include = True
    assert len(config.tabs.children) == 3
    parameters = model.get_model_state()
    assert "bands" in parameters


def test_reminder_info():
    """Dynamic add/remove the reminder text based on the workchain settings."""
    model = ConfigurationModel()
    config = ConfigureQeAppWorkChainStep(model=model)
    config.render()
    bands_info = config.property_children[0].children[1]
    assert bands_info.value == ""
    bands_model = model.get_model("bands")
    bands_model.include = True
    assert bands_info.value == "Customize bands settings in <b>Step 2.2</b> if needed"
    bands_model.include = False
    assert bands_info.value == ""
