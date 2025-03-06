import tempfile
from pathlib import Path

from aiidalab_qe.app.utils.plugin_manager import PluginManager

# mock the content of the YAML file
yaml_content = """
    my-plugin:
      title: My Test Plugin
      description: A test plugin
      author: John Doe
      pip: my-plugin
      status: beta
    another-plugin:
      title: Another Plugin
      description: Another test plugin
      author: Jane Doe
      pip: another-plugin
      status: experimental
    """


def test_plugin_manager_local_config_file():
    """
    Test that PluginManager loads the YAML file and builds the UI (Accordion) properly.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        config_file = tmp_path / "test_plugins.yaml"
        config_file.write_text(yaml_content)

        # Instantiate the manager
        manager = PluginManager(config_source=str(config_file))
        manager._build_ui()
        assert len(manager.accordion.children) == 2, "Should build 2 accordion panels."

        # Check that the titles are set correctly
        assert "My Test Plugin" in manager.accordion.get_title(
            0
        ), "Check initial title for not-installed plugin"
        assert "Another Plugin" in manager.accordion.get_title(
            1
        ), "Check second plugin title"


def test_plugin_manager_default_config_file():
    """
    Test that PluginManager loads the YAML file from default source (GitHub repo).
    """
    manager = PluginManager()
    manager._build_ui()
    assert (
        len(manager.accordion.children) > 0
    ), "Should build at least one accordion panels."
