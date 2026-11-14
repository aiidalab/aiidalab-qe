from aiidalab_qe.plugins import QeAppPlugin
from aiidalab_qe.plugins.utils import PluginCatalog


class ExamplePlugin(QeAppPlugin):
    @property
    def outline(self) -> str:
        return "Example plugin"

    def configuration(self):
        return {"panel": object, "model": object}


def test_plugin_catalog_loads_requested_features():
    plugin = ExamplePlugin()
    catalog = PluginCatalog({"example": plugin})

    assert catalog.plugins["example"] is plugin
    assert catalog.get_features("outline") == {"example": "Example plugin"}
    assert catalog.get_features("configuration") == {
        "example": {"panel": object, "model": object}
    }


def test_optional_features_default_to_absent():
    class OutlineOnlyPlugin(QeAppPlugin):
        @property
        def outline(self) -> str:
            return "Outline only"

    plugin = OutlineOnlyPlugin()

    for feature in ("configuration", "resources", "workchain", "examples"):
        assert plugin.load_feature(feature) is None


def test_plugin_catalog_filters_before_feature_load():
    class UnavailablePlugin(ExamplePlugin):
        def configuration(self):
            raise AssertionError("filtered plugin feature must not load")

    catalog = PluginCatalog(
        {"available": ExamplePlugin(), "unavailable": UnavailablePlugin()}
    )

    features = catalog.get_features(
        "configuration",
        plugin_filter=lambda name, _plugin: name == "available",
    )

    assert list(features) == ["available"]


def test_plugin_importer_and_editor_features_can_return_multiple_components():
    class MultiComponentPlugin(ExamplePlugin):
        def importers(self):
            return (str, bytes)

        def editors(self):
            return (int, float)

    catalog = PluginCatalog({"multi": MultiComponentPlugin()})

    assert catalog.get_features("importers") == {"multi": (str, bytes)}
    assert catalog.get_features("editors") == {"multi": (int, float)}


def test_plugin_catalog_reports_feature_load_failure(capsys):
    class BrokenPlugin(ExamplePlugin):
        def configuration(self):
            raise RuntimeError("configuration import failed")

    catalog = PluginCatalog({"broken": BrokenPlugin()})

    assert catalog.get_features("configuration") == {}
    assert "Failed to load plugin feature configuration for broken." in (
        capsys.readouterr().out
    )


def test_plugin_registration_accepts_plugin_class(monkeypatch):
    from aiidalab_qe.plugins import utils
    from aiidalab_qe.plugins.utils import test_plugin_functionality as validate_plugin

    class EntryPoint:
        def load(self):
            return ExamplePlugin

    monkeypatch.setattr(
        utils,
        "get_plugin_entry_points_for_package",
        lambda _package: [EntryPoint()],
    )

    assert validate_plugin("example-plugin") == (True, "")


def test_plugin_registration_rejects_non_plugin_class(monkeypatch):
    from aiidalab_qe.plugins import utils
    from aiidalab_qe.plugins.utils import test_plugin_functionality as validate_plugin

    class NotAPlugin:
        pass

    class EntryPoint:
        def load(self):
            return NotAPlugin

    monkeypatch.setattr(
        utils,
        "get_plugin_entry_points_for_package",
        lambda _package: [EntryPoint()],
    )

    success, message = validate_plugin("example-plugin")
    assert not success
    assert "QeAppPlugin subclass" in message
