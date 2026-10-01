import tempfile
from io import StringIO
from pathlib import Path

import ipywidgets as ipw
import pytest

from aiidalab_qe.app.utils import plugin_manager
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
        assert "My Test Plugin" in manager.accordion.get_title(0), (
            "Check initial title for not-installed plugin"
        )
        assert "Another Plugin" in manager.accordion.get_title(1), (
            "Check second plugin title"
        )


def test_outdated_plugin_actions(monkeypatch):
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.0", False, None),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)

    with tempfile.TemporaryDirectory() as tmp_dir:
        config_file = Path(tmp_dir) / "test_plugins.yaml"
        config_file.write_text(yaml_content)
        manager = PluginManager(config_source=str(config_file))
        manager._build_ui()

    buttons = manager.accordion.children[0].children[2].children
    install_button, update_button, _post_install_button, remove_button, clear_button = (
        buttons
    )
    assert install_button.disabled
    assert not update_button.disabled
    assert not remove_button.disabled
    assert clear_button.description == "Clear output"
    warning = manager.accordion.children[0].children[1].value
    assert 'class="alert alert-danger"' in warning
    assert "Installed version 1.0" in warning


def test_install_button_stays_disabled_after_update(monkeypatch):
    class Process:
        returncode = 0
        stdout = StringIO()

        def poll(self):
            return self.returncode

    monkeypatch.setattr(
        plugin_manager.subprocess, "Popen", lambda *_args, **_kwargs: Process()
    )
    install_button = ipw.Button()
    remove_button = ipw.Button(disabled=True)

    result = plugin_manager.execute_command_with_output(
        ["pip", "install", "--upgrade", "my-plugin"],
        ipw.HTML(),
        install_button,
        remove_button,
        action="update",
    )

    assert result
    assert install_button.disabled
    assert not remove_button.disabled


def test_clear_output_button_clears_and_hides_logs(monkeypatch):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {"my-plugin": {"title": "My Test Plugin", "pip": "my-plugin"}},
    )
    monkeypatch.setattr(
        plugin_manager, "get_plugin_version_info", lambda *_args: (None, True, None)
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)

    manager = PluginManager()
    manager._build_ui()
    row = manager.accordion.children[0].children
    message_container, output_container = row[3:5]
    clear_button = row[2].children[-1]
    assert clear_button.disabled

    message_container.value = "Status"
    assert not clear_button.disabled

    output_container.value = "Command output"
    message_container.layout.display = "block"
    output_container.layout.display = "block"

    clear_button.click()

    assert message_container.value == ""
    assert output_container.value == ""
    assert clear_button.disabled
    assert message_container.layout.display == "none"
    assert output_container.layout.display == "none"


@pytest.mark.parametrize(
    ("plugin_compatible", "app_compatible", "expected_icon"),
    [(True, True, "✅"), (False, True, "⚠️"), (True, False, "⚠️")],
)
def test_installed_plugin_status_icon(
    monkeypatch, plugin_compatible, app_compatible, expected_icon
):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "pip": "my-plugin>=1",
            }
        },
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.0", plugin_compatible, None),
    )
    monkeypatch.setattr(
        plugin_manager, "is_version_compatible", lambda *_args: app_compatible
    )

    manager = PluginManager()
    manager._build_ui()

    assert manager.accordion.get_title(0).endswith(expected_icon)


def test_plugin_manager_displays_accordion_without_refresh_button(monkeypatch):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {"my-plugin": {"title": "My Test Plugin", "pip": "my-plugin"}},
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: (None, True, None),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)
    displayed = []
    monkeypatch.setattr(plugin_manager, "display", displayed.append)

    manager = PluginManager()
    manager.display_ui()

    assert displayed == [manager.accordion]
    assert not hasattr(manager, "refresh_button")


def test_plugin_manager_default_config_file():
    """
    Test that PluginManager loads the YAML file from default source (GitHub repo).
    """
    manager = PluginManager()
    manager._build_ui()
    assert len(manager.accordion.children) > 0, (
        "Should build at least one accordion panels."
    )


@pytest.mark.parametrize(
    ("installed_version", "requirement", "expected"),
    [
        ("1.2.8", "my-plugin>=1.2.9", False),
        ("1.2.9", "my-plugin>=1.2.9", True),
        ("1.3.0", "my-plugin>=1.2.9", True),
        ("1.0", "my-plugin", True),
    ],
)
def test_get_plugin_version_info(monkeypatch, installed_version, requirement, expected):
    monkeypatch.setattr(
        plugin_manager.metadata, "version", lambda _package: installed_version
    )

    version, compatible, error = plugin_manager.get_plugin_version_info(
        "my-plugin", requirement
    )

    assert version == installed_version
    assert compatible is expected
    assert error is None


def test_get_plugin_version_info_not_installed(monkeypatch):
    def missing_package(_package):
        raise plugin_manager.metadata.PackageNotFoundError

    monkeypatch.setattr(plugin_manager.metadata, "version", missing_package)

    assert plugin_manager.get_plugin_version_info("my-plugin", "my-plugin>=1") == (
        None,
        True,
        None,
    )


def test_get_plugin_version_info_invalid_requirement():
    version, compatible, error = plugin_manager.get_plugin_version_info(
        "my-plugin", "not a requirement"
    )

    assert version is None
    assert not compatible
    assert error


def test_update_package_clears_warning_after_minimum_is_met(monkeypatch):
    commands = []
    monkeypatch.setattr(
        plugin_manager,
        "execute_command_with_output",
        lambda command, *_args, **_kwargs: commands.append(command) or True,
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    install_button = ipw.Button()
    update_button = ipw.Button()
    remove_button = ipw.Button()
    warning = ipw.HTML(value="Outdated")
    accordion = ipw.Accordion(children=[ipw.VBox()])
    accordion.set_title(0, "My Test Plugin ⚠️")

    plugin_manager.update_package(
        "my-plugin",
        "my-plugin>=1.2.9",
        "",
        ipw.HTML(),
        ipw.HTML(),
        install_button,
        update_button,
        remove_button,
        warning,
        accordion,
        0,
    )

    assert commands == [["pip", "install", "--upgrade", "my-plugin>=1.2.9", "--user"]]
    assert warning.value == ""
    assert update_button.disabled
    assert accordion.get_title(0).endswith("✅")


def test_update_package_keeps_warning_if_minimum_is_not_met(monkeypatch):
    monkeypatch.setattr(
        plugin_manager,
        "execute_command_with_output",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.2.8", False, None),
    )
    warning = ipw.HTML(value="Outdated")
    message = ipw.HTML()
    update_button = ipw.Button()

    plugin_manager.update_package(
        "my-plugin",
        "my-plugin>=1.2.9",
        "",
        ipw.HTML(),
        message,
        ipw.Button(),
        update_button,
        ipw.Button(),
        warning,
    )

    assert warning.value == "Outdated"
    assert not update_button.disabled
    assert "does not meet the required version" in message.value


def test_run_post_install_runs_only_configured_command(monkeypatch):
    calls = []

    class Result:
        returncode = 0
        stdout = "setup complete"
        stderr = ""

    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **kwargs: calls.append(command) or Result(),
    )
    output = ipw.HTML()
    message = ipw.HTML()

    plugin_manager.run_post_install("my-plugin", "setup", output, message)

    assert calls == [[plugin_manager.sys.executable, "-m", "my_plugin", "setup"]]
    assert "background-color: #3B3B3B" in output.value
    assert "color: #FFFFFF" in output.value
    assert "setup complete" in output.value
    assert "Post-install completed" in message.value


def test_post_install_only_button_is_enabled_when_package_is_not_installed(
    monkeypatch,
):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "pip": "my-plugin",
                "post_install": "setup",
            }
        },
    )
    monkeypatch.setattr(
        plugin_manager, "get_plugin_version_info", lambda *_args: (None, True, None)
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)

    manager = PluginManager()
    manager._build_ui()

    post_install_button = manager.accordion.children[0].children[2].children[2]
    assert not post_install_button.disabled
