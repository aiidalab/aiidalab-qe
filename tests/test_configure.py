from aiidalab_qe.app.configuration import ConfigurationStep, ConfigurationStepModel
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION


def test_get_configuration_parameters():
    model = ConfigurationStepModel()
    _ = ConfigurationStep(model=model)
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
    model = ConfigurationStepModel()
    _ = ConfigurationStep(model=model)
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
    model = ConfigurationStepModel()
    config = ConfigurationStep(model=model)
    config.render()
    assert len(config.tabs.children) == 2
    parameters = model.get_model_state()
    assert "bands" not in parameters
    model.get_model("bands").include = True
    assert len(config.tabs.children) == 3
    parameters = model.get_model_state()
    assert "bands" in parameters


def test_fetching_available_properties():
    import os

    current_file = os.path.abspath(__file__)
    plugin_file = os.path.join(os.path.dirname(current_file), "../plugins.yaml")
    model = ConfigurationStepModel()
    config = ConfigurationStep(model=model)
    config._fetch_available_properties(str(plugin_file))
    assert len(config.available_properties_list) > 0
    assert model.available_properties_fetched
    assert all("<li>" not in title for title in config.available_properties_list)


def test_available_properties_render_with_single_list_item(monkeypatch):
    from aiidalab_qe.app.configuration import step as configuration_step

    class Registry:
        data = {
            "my-plugin": {
                "title": "My plugin",
                "pip": "my-plugin",
                "category": "calculation",
            }
        }

        def __init__(self, _source):
            pass

    monkeypatch.setattr(configuration_step, "PluginManager", Registry)
    monkeypatch.setattr(
        configuration_step,
        "get_plugin_version_info",
        lambda *_args: (None, True, None),
    )
    monkeypatch.setattr(
        configuration_step, "is_version_compatible", lambda *_args: True
    )

    config = ConfigurationStep(model=ConfigurationStepModel())
    config.render()

    assert config.available_properties_list == ["My plugin"]
    assert config.available_properties.value.count("<li>") == 1
    assert "<li></li>" not in config.available_properties.value


def test_incompatible_plugin_is_listed_and_not_restored(monkeypatch):
    from types import SimpleNamespace

    from aiidalab_qe.app.configuration import step as configuration_step

    class Registry:
        data = {
            "my-plugin": {
                "title": "My plugin",
                "package": "my-plugin",
                "pip": "my-plugin>=2.0",
                "category": "calculation",
            }
        }

        def __init__(self, _source):
            pass

    distribution = SimpleNamespace(
        metadata={"Name": "my-plugin"},
        entry_points=[SimpleNamespace(group="aiidalab_qe.properties", name="bands")],
    )
    monkeypatch.setattr(configuration_step, "PluginManager", Registry)
    monkeypatch.setattr(
        configuration_step,
        "distributions",
        lambda: [distribution],
    )
    monkeypatch.setattr(
        configuration_step,
        "get_plugin_version_info",
        lambda *_args: ("1.0", False, None),
    )
    monkeypatch.setattr(
        configuration_step, "is_version_compatible", lambda *_args: True
    )

    model = ConfigurationStepModel()
    config = ConfigurationStep(model=model)
    config.render()

    assert config.incompatible_properties_list == ["My plugin"]
    assert "Electronic band structure" not in [
        row.children[0].title for row in config.installed_properties_list
    ]

    model.set_model_state(
        {"workchain": {"properties": ["bands"], "relax_type": "none"}}
    )
    assert not model.get_model("bands").include
    assert "bands" not in model._get_properties()
