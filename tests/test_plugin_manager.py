from aiidalab_qe.app.utils.plugin_manager import PluginManager


def test_plugin_manager_init():
    """
    Test that PluginManager loads the YAML file and builds the UI (Accordion) properly.
    """
    # Instantiate the manager
    manager = PluginManager()
    assert len(manager.accordion.children) == 2, "Should build 2 accordion panels."
    # Check that the titles are set correctly
    assert "My Test Plugin" in manager.accordion.get_title(
        0
    ), "Check initial title for not-installed plugin"
    assert "Another Plugin" in manager.accordion.get_title(
        1
    ), "Check second plugin title"
