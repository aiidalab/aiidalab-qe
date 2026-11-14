def test_entry_point():
    from aiidalab_qe.plugins.utils import get_plugins, load_plugin_catalog

    load_plugin_catalog.cache_clear()
    entries = get_plugins()
    assert "bands" in entries
    assert "pdos" in entries
    assert entries["bands"].workchain() is not None

    entries_list = list(entries.keys())
    prioritized_entries = ["electronic_structure", "bands", "pdos"]
    for i, prioritized_entry in enumerate(prioritized_entries):
        assert entries_list.index(prioritized_entry) == i, (
            f"Entry point {prioritized_entry} is not in the expected position."
        )


def test_entry_point_instantiates_qe_app_plugin_lazily(monkeypatch, request):
    from aiidalab_qe.plugins import QeAppPlugin
    from aiidalab_qe.plugins.utils import (
        get_plugin_feature,
        get_plugins,
        load_plugin_catalog,
    )

    loaded_features = []

    class ExamplePlugin(QeAppPlugin):
        @property
        def outline(self):
            return "Example plugin"

        def configuration(self):
            loaded_features.append("configuration")
            return {"panel": object, "model": object}

    class EntryPoint:
        name = "example"

        def load(self):
            return ExamplePlugin

    from aiidalab_qe.plugins import utils

    load_plugin_catalog.cache_clear()
    request.addfinalizer(load_plugin_catalog.cache_clear)
    monkeypatch.setattr(utils, "entry_points", lambda **_: [EntryPoint()])

    plugin = get_plugins()["example"]
    assert isinstance(plugin, ExamplePlugin)
    assert get_plugins()["example"] is plugin
    assert loaded_features == []
    assert plugin.outline == "Example plugin"
    assert get_plugin_feature(plugin, "configuration") == {
        "panel": object,
        "model": object,
    }
    assert loaded_features == ["configuration"]
