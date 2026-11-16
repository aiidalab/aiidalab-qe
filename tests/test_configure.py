from aiidalab_qe.app.configuration import ConfigurationStep, ConfigurationStepModel
from aiidalab_qe.plugins import registry as plugin_registry
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION


def _load_test_registry(monkeypatch, plugins):
    monkeypatch.setattr(plugin_registry, "load_plugin_config", lambda: plugins)
    return plugin_registry.PluginRegistry.load()


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
    assert config.plugin_validation_message_container.layout.margin == "0"
    assert config.plugin_validation_message_container.layout.display == "none"
    assert config.plugin_validation_message_container.children == ()
    assert (
        config.plugin_sections.children[2] is config.plugin_validation_message_container
    )
    parameters = model.get_model_state()
    assert "bands" not in parameters
    model.get_model("bands").include = True
    assert len(config.tabs.children) == 3
    parameters = model.get_model_state()
    assert "bands" in parameters


def test_plugin_validation_message_is_shown_and_hidden(monkeypatch):
    from aiidalab_widgets_base import LoadingWidget

    config = ConfigurationStep(model=ConfigurationStepModel())
    monkeypatch.setattr(config, "_schedule_configuration_plugin_loading", lambda: None)
    config.render()

    config._set_plugin_validation_message("Checking plugin compatibility")
    assert isinstance(config.plugin_validation_loading_widget, LoadingWidget)
    assert (
        config.plugin_validation_loading_widget.message.value
        == "Checking plugin compatibility"
    )
    assert config.plugin_validation_message_container.children == (
        config.plugin_validation_loading_widget,
    )
    assert config.plugin_validation_message_container.layout.display == ""

    config._set_plugin_validation_message(None)
    assert config.plugin_validation_message_container.layout.display == "none"
    assert config.plugin_validation_message_container.children == ()

    config._set_plugin_validation_message("Plugin validation failed", error=True)
    assert "Plugin validation failed" in (
        config.plugin_validation_message_container.children[0].value
    )


def test_fetching_available_plugins(monkeypatch):
    from types import SimpleNamespace

    model = ConfigurationStepModel()
    config = ConfigurationStep(model=model)
    monkeypatch.setattr(config, "_schedule_plugin_status_check", lambda: None)
    config.render()
    statuses = {
        name: SimpleNamespace(is_installed=False, is_incompatible=False)
        for name in config._plugin_registry.plugins
    }
    config._update_plugin_lists(statuses)
    assert len(config.available_plugins_list) > 0
    assert all("<li>" not in title for title in config.available_plugins_list)


def test_available_plugins_render_with_single_list_item(monkeypatch):
    from aiidalab_qe.app.configuration import step as configuration_step

    registry = _load_test_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My plugin",
                "description": "Plugin description",
                "pip": "my-plugin",
                "category": "calculation",
            }
        },
    )
    monkeypatch.setattr(configuration_step, "get_plugin_registry", lambda: registry)
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (None, True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

    config = ConfigurationStep(model=ConfigurationStepModel())
    config.render()

    assert config.available_plugins_list == ["My plugin"]
    available_plugins_html = config.available_plugins.children[0].value
    assert available_plugins_html.count("<li>") == 1
    assert "<li></li>" not in available_plugins_html
    assert config.available_plugins_subtitle.value == (
        "<em>The following plugins are available in the plugin store:</em>"
    )


def test_available_plugins_subtitle_when_none_are_available():
    import ipywidgets as ipw

    config = ConfigurationStep(model=ConfigurationStepModel())
    config.available_plugins = ipw.VBox()
    config.available_plugins_subtitle = ipw.HTML()
    config.available_plugins_list = []

    config._render_available_plugins_list()

    assert config.available_plugins_subtitle.value == (
        "<em>All available plugins are installed</em>"
    )


def test_incompatible_plugin_is_listed_and_not_restored(monkeypatch):
    from types import SimpleNamespace

    from aiidalab_qe.app.configuration import step as configuration_step

    distribution = SimpleNamespace(
        metadata={"Name": "my-plugin"},
        entry_points=[SimpleNamespace(group="aiidalab_qe.plugins", name="bands")],
    )
    registry = _load_test_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My plugin",
                "description": "Plugin description",
                "pip": "my-plugin>=2.0",
                "category": "calculation",
            }
        },
    )
    monkeypatch.setattr(configuration_step, "get_plugin_registry", lambda: registry)
    monkeypatch.setattr(
        plugin_registry.metadata, "distributions", lambda: [distribution]
    )
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.0", False, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

    model = ConfigurationStepModel()
    config = ConfigurationStep(model=model)
    config.render()

    assert config.incompatible_plugins_list == ["My plugin"]
    assert config.installed_plugin_rows["bands"] not in config.installed_plugins_list

    model.set_model_state(
        {"workchain": {"properties": ["bands"], "relax_type": "none"}}
    )
    assert not model.get_model("bands").include
    assert "bands" not in model._get_properties()


def test_activation_failed_plugin_is_not_ready_for_new_calculations(monkeypatch):
    from types import SimpleNamespace

    from aiidalab_qe.app.configuration import step as configuration_step

    distribution = SimpleNamespace(
        metadata={"Name": "my-plugin"},
        entry_points=[SimpleNamespace(group="aiidalab_qe.plugins", name="bands")],
    )
    registry = _load_test_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My plugin",
                "description": "Plugin description",
                "pip": "my-plugin>=2.0",
                "category": "calculation",
            }
        },
    )
    monkeypatch.setattr(configuration_step, "get_plugin_registry", lambda: registry)
    monkeypatch.setattr(
        plugin_registry.metadata, "distributions", lambda: [distribution]
    )
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("2.0", True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_: True)
    monkeypatch.setattr(
        plugin_registry,
        "poll_plugin_registration",
        lambda _package: "Plugin validation failed",
    )

    config = ConfigurationStep(model=ConfigurationStepModel())
    config.render()

    assert config.incompatible_plugins_list == ["My plugin"]
    assert config.installed_plugin_rows["bands"] not in config.installed_plugins_list


def test_process_loaded_plugin_is_ready_and_restored(monkeypatch):
    from types import SimpleNamespace

    from aiidalab_qe.app.configuration import step as configuration_step

    distribution = SimpleNamespace(
        metadata={"Name": "my-plugin"},
        entry_points=[SimpleNamespace(group="aiidalab_qe.plugins", name="bands")],
    )
    registry = _load_test_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My plugin",
                "description": "Plugin description",
                "pip": "my-plugin>=2.0",
                "category": "calculation",
            }
        },
    )
    monkeypatch.setattr(configuration_step, "get_plugin_registry", lambda: registry)
    monkeypatch.setattr(
        plugin_registry.metadata, "distributions", lambda: [distribution]
    )
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (_ for _ in ()).throw(
            AssertionError("process-backed Step 2 must skip compatibility checks")
        ),
    )
    monkeypatch.setattr(
        plugin_registry,
        "is_version_compatible",
        lambda *_args: (_ for _ in ()).throw(
            AssertionError("process-backed Step 2 must skip compatibility checks")
        ),
    )
    monkeypatch.setattr(plugin_registry, "is_plugin_installed", lambda *_: True)
    model = ConfigurationStepModel()
    model.loaded_from_process = True
    config = ConfigurationStep(model=model)
    config.render()
    model.set_model_state(
        {"workchain": {"properties": ["bands"], "relax_type": "none"}}
    )

    assert config.incompatible_plugins_list == []
    assert config.installed_plugin_rows["bands"] in config.installed_plugins_list
    assert model.get_model("bands").include
    assert "bands" in model._get_properties()
