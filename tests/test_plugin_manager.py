import tempfile
from io import StringIO
from pathlib import Path

import ipywidgets as ipw
import pytest

from aiidalab_qe.app.utils import plugin_manager
from aiidalab_qe.app.utils.plugin_manager import (
    PluginManager,
    QeAppPlugin,
    QeAppPluginData,
)
from aiidalab_qe.plugins import state as plugin_state


@pytest.fixture(autouse=True)
def isolate_plugin_activation_state(monkeypatch, tmp_path):
    monkeypatch.setattr(
        plugin_state,
        "ACTIVATION_STATE_PATH",
        tmp_path / "plugin-activation.json",
    )


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


def test_checked_in_plugin_registry_is_valid():
    registry = plugin_manager.yaml.safe_load(
        plugin_manager.DEFAULT_PLUGIN_CONFIG_SOURCE.read_text()
    )

    plugins = [
        QeAppPluginData.from_mapping(name, data) for name, data in registry.items()
    ]

    assert plugins


@pytest.mark.parametrize(
    ("data", "message"),
    [
        ({"description": "missing title", "pip": "my-plugin"}, "title"),
        (
            {"title": "My plugin", "description": "Test", "pip": "not a requirement"},
            "invalid 'pip' requirement",
        ),
        (
            {
                "title": "My plugin",
                "description": "Test",
                "github": "https://example.test/plugin",
                "requires_aiidalab_qe": "invalid",
            },
            "invalid 'requires_aiidalab_qe' specifier",
        ),
    ],
)
def test_plugin_registry_entry_validation(data, message):
    with pytest.raises(ValueError, match=message):
        QeAppPluginData.from_mapping("my-plugin", data)


def test_plugin_registry_preserves_unknown_metadata():
    plugin = QeAppPluginData.from_mapping(
        "my-plugin",
        {
            "title": "My plugin",
            "description": "Test",
            "pip": "my-plugin",
            "future_field": "kept",
        },
    )

    assert plugin.as_mapping()["future_field"] == "kept"


def test_plugin_data_and_widget_keep_metadata_and_widgets_separate():
    data = QeAppPluginData.from_mapping(
        "my-plugin",
        {
            "title": "My plugin",
            "description": "Test",
            "pip": "my-plugin",
        },
    )

    plugin = QeAppPlugin(data, plugin_name="my-plugin")

    assert plugin.data is data
    assert isinstance(plugin, ipw.VBox)
    assert not hasattr(data, "version_warning")
    assert not hasattr(data, "install_button")
    assert isinstance(plugin.version_warning, ipw.HTML)
    assert isinstance(plugin.install_button, ipw.Button)


def test_empty_registry_is_normalized_to_mapping(tmp_path):
    config_file = tmp_path / "plugins.yaml"
    config_file.write_text("")

    manager = PluginManager(config_source=str(config_file))

    assert manager.data == {}
    assert manager.config_error is None


def test_invalid_registry_entry_is_reported_in_store(tmp_path, monkeypatch):
    config_file = tmp_path / "plugins.yaml"
    config_file.write_text(
        "my-plugin:\n  title: My plugin\n  description: Test\n  pip: 'not a requirement'\n"
    )
    displayed = []
    monkeypatch.setattr(plugin_manager, "display", displayed.append)

    manager = PluginManager(config_source=str(config_file))
    manager.display_ui()

    assert "invalid 'pip' requirement" in manager.config_error
    assert "invalid &#x27;pip&#x27; requirement" in displayed[0].value


def test_remote_registry_request_error_is_reported(monkeypatch):
    def fail_request(*_args, **_kwargs):
        raise plugin_manager.requests.ConnectionError("offline")

    monkeypatch.setattr(plugin_manager.requests, "get", fail_request)

    manager = PluginManager(config_source="https://example.test/plugins.yaml")

    assert "Could not fetch plugin registry: offline" in manager.config_error


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
    (
        install_button,
        _post_install_button,
        update_button,
        remove_button,
        _retry_activation_button,
        clear_button,
    ) = buttons
    assert install_button.disabled
    assert not update_button.disabled
    assert not remove_button.disabled
    assert clear_button.description == "Clear output"
    warning = manager.accordion.children[0].children[1].value
    assert 'class="alert alert-danger"' in warning
    assert "Installed version 1.0" in warning


def test_execute_command_streams_without_changing_widget_policy(monkeypatch):
    class Process:
        returncode = 0
        stdout = StringIO()

        def wait(self):
            return self.returncode

    monkeypatch.setattr(
        plugin_manager.subprocess, "Popen", lambda *_args, **_kwargs: Process()
    )
    row = QeAppPlugin(
        QeAppPluginData.from_mapping(
            "my-plugin",
            {
                "title": "My plugin",
                "description": "Test",
                "pip": "my-plugin",
            },
        ),
        plugin_name="my-plugin",
    )
    row.output_container = ipw.HTML()

    result = row._execute_command(["pip", "install", "--upgrade", "my-plugin"])

    assert result
    assert row.output_container.value == ""


def test_clear_output_button_clears_and_hides_logs(monkeypatch):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin",
            }
        },
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
                "description": "A test plugin",
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


def test_plugin_manager_displays_persistent_ui_with_refresh_button(monkeypatch):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin",
            }
        },
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
    first_row = manager.accordion.children[0]
    first_row.message_container.value = "Keep this log"
    manager.display_ui()

    assert displayed == [manager.ui, manager.ui]
    assert manager.ui.children == (manager.refresh_button, manager.accordion)
    assert manager.accordion.children[0] is first_row
    assert first_row.message_container.value == "Keep this log"


def test_manager_refresh_updates_existing_rows_without_clearing_logs(monkeypatch):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            name: {
                "title": name,
                "description": "A test plugin",
                "pip": name,
            }
            for name in ("first-plugin", "second-plugin")
        },
    )
    version_info = {
        "first-plugin": (None, True, None),
        "second-plugin": (None, True, None),
    }
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda name, *_args: version_info[name],
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)

    manager = PluginManager()
    manager._build_ui()
    first_row, second_row = manager.accordion.children
    first_row.message_container.value = "Previous update log"
    first_row.output_container.value = "Previous command output"
    plugin_state.set_activation_failure("first-plugin", "External activation failure")
    version_info["first-plugin"] = ("1.0", False, None)
    version_info["second-plugin"] = ("2.0", True, None)

    manager.refresh_button.click()

    assert manager.accordion.children == (first_row, second_row)
    assert first_row.installed_version == "1.0"
    assert not first_row.plugin_compatible
    assert first_row.activation_error == "External activation failure"
    assert "Installed version 1.0" in first_row.version_warning.value
    assert "External activation failure" in first_row.version_warning.value
    assert first_row.version_warning_text == first_row.version_warning.value
    assert manager.accordion.get_title(0).endswith("⚠️")
    assert not first_row.update_button.disabled
    assert first_row.message_container.value == "Previous update log"
    assert first_row.output_container.value == "Previous command output"
    assert second_row.installed_version == "2.0"


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
        ("1.3.0", "my-plugin>=1.2.9,<1.3.0", False),
        ("1.2.9", "my-plugin>=1.2.9,<1.3.0", True),
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


def test_plugin_activation_failure_persists_and_clears(monkeypatch, tmp_path):
    from aiidalab_qe.plugins import state

    state_path = tmp_path / "plugin-activation.json"
    monkeypatch.setattr(state, "ACTIVATION_STATE_PATH", state_path)

    state.set_activation_failure("my-plugin", "Plugin validation failed")

    assert state.get_activation_failure("my-plugin") == "Plugin validation failed"
    assert state_path.exists()

    state.set_activation_failure("another-plugin", "Daemon restart failed")
    state.clear_activation_failure("my-plugin")

    assert state.get_activation_failure("my-plugin") is None
    assert state.get_activation_failure("another-plugin") == "Daemon restart failed"

    state.clear_activation_failure("another-plugin")
    assert not state_path.exists()


@pytest.mark.parametrize("contents", ["not valid json", "[]"])
def test_invalid_plugin_activation_state_is_logged_and_recovered(
    monkeypatch, tmp_path, caplog, contents
):
    state_path = tmp_path / "plugin-activation.json"
    state_path.write_text(contents, encoding="utf-8")
    monkeypatch.setattr(plugin_state, "ACTIVATION_STATE_PATH", state_path)

    assert plugin_state.get_activation_failure("my-plugin") is None
    plugin_state.clear_activation_failure("my-plugin")
    plugin_state.set_activation_failure("my-plugin", "Validation failed")

    assert plugin_state.get_activation_failure("my-plugin") == "Validation failed"
    assert "plugin activation state" in caplog.text.lower()


def test_update_package_clears_warning_after_minimum_is_met(monkeypatch):
    commands = []

    def execute(_row, command, *_args, **_kwargs):
        commands.append(command)
        installed_version["value"] = "1.2.9"
        return True

    installed_version = {"value": "1.2.8"}
    monkeypatch.setattr(QeAppPlugin, "_execute_command", execute)
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: (
            installed_version["value"],
            installed_version["value"] == "1.2.9",
            None,
        ),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)
    daemon_commands = []
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: (
            daemon_commands.append(command)
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin>=1.2.9",
                "post_install": "setup",
            }
        },
    )
    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]
    plugin._on_update(None)

    assert commands == [
        [
            plugin_manager.sys.executable,
            "-m",
            "pip",
            "install",
            "--upgrade",
            "my-plugin>=1.2.9",
            "--user",
        ]
    ]
    assert plugin.version_warning.value == ""
    assert plugin.update_button.disabled
    assert manager.accordion.get_title(0).endswith("✅")
    assert daemon_commands == [
        [plugin_manager.sys.executable, "-m", "my_plugin", "setup"],
        [
            plugin_manager.sys.executable,
            "-m",
            "aiidalab_qe",
            "test-plugin",
            "my-plugin",
        ],
        ["verdi", "daemon", "restart"],
    ]


def test_update_package_resolves_version_above_upper_bound(monkeypatch):
    commands = []
    installed_version = {"value": "1.3.0"}
    requirement = "my-plugin>=1.2.9,<1.3.0"

    def execute(_row, command, *_args, **_kwargs):
        commands.append(command)
        installed_version["value"] = "1.2.9"
        return True

    monkeypatch.setattr(QeAppPlugin, "_execute_command", execute)
    monkeypatch.setattr(
        plugin_manager.metadata,
        "version",
        lambda _name: installed_version["value"],
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: (
            installed_version["value"],
            installed_version["value"] == "1.2.9",
            None,
        ),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)
    subprocess_commands = []
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: (
            subprocess_commands.append(command)
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": requirement,
            }
        },
    )
    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]
    plugin._on_update(None)

    assert commands == [
        [
            plugin_manager.sys.executable,
            "-m",
            "pip",
            "install",
            "--upgrade",
            requirement,
            "--user",
        ]
    ]
    assert plugin.version_warning.value == ""
    assert plugin.update_button.disabled
    assert subprocess_commands == [
        [
            plugin_manager.sys.executable,
            "-m",
            "aiidalab_qe",
            "test-plugin",
            "my-plugin",
        ],
        ["verdi", "daemon", "restart"],
    ]


def test_update_does_not_restart_when_plugin_test_fails(monkeypatch):
    installed_version = {"value": "1.2.8"}
    commands = []

    def execute(_row, command, *_args, **_kwargs):
        commands.append(command)
        installed_version["value"] = "1.2.9"
        return True

    monkeypatch.setattr(QeAppPlugin, "_execute_command", execute)
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: (
            installed_version["value"],
            installed_version["value"] == "1.2.9",
            None,
        ),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)

    subprocess_commands = []

    def fail_plugin_test(command, **_kwargs):
        subprocess_commands.append(command)
        return plugin_manager.subprocess.CompletedProcess(command, 1, "", "failed")

    monkeypatch.setattr(plugin_manager.subprocess, "run", fail_plugin_test)
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin>=1.2.9",
            }
        },
    )
    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]

    plugin._on_update(None)

    assert len(commands) == 1
    assert subprocess_commands == [
        [
            plugin_manager.sys.executable,
            "-m",
            "aiidalab_qe",
            "test-plugin",
            "my-plugin",
        ]
    ]
    assert "did not pass" in plugin.message_container.value
    assert "plugin remains installed" in plugin.message_container.value
    assert "daemon was not restarted" in plugin.message_container.value
    assert "Updated my-plugin" not in plugin.message_container.value
    assert plugin.activation_error
    assert plugin.retry_activation_button.layout.display == ""
    assert plugin.update_button.disabled is False
    assert plugin_state.get_activation_failure("my-plugin") == plugin.activation_error


def test_retry_activation_clears_failure_and_restarts_daemon(monkeypatch):
    plugin_state.set_activation_failure("my-plugin", "Previous validation failed")
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)
    calls = []

    def run(command, **_kwargs):
        calls.append(command)
        return plugin_manager.subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(plugin_manager.subprocess, "run", run)
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin>=1.2.9",
            }
        },
    )
    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]

    assert plugin.retry_activation_button.layout.display == ""
    plugin._on_retry_activation(None)

    assert calls == [
        [
            plugin_manager.sys.executable,
            "-m",
            "aiidalab_qe",
            "test-plugin",
            "my-plugin",
        ],
        ["verdi", "daemon", "restart"],
    ]
    assert plugin_state.get_activation_failure("my-plugin") is None
    assert plugin.retry_activation_button.layout.display == "none"
    assert manager.accordion.get_title(0).endswith("✅")


@pytest.mark.parametrize("failure_stage", ["post_install", "daemon_restart"])
def test_update_failure_keeps_activation_warning(monkeypatch, failure_stage):
    monkeypatch.setattr(
        QeAppPlugin,
        "_execute_command",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)
    commands = []

    def run(command, **_kwargs):
        commands.append(command)
        if command[1:3] == ["-m", "my_plugin"]:
            return plugin_manager.subprocess.CompletedProcess(
                command,
                1 if failure_stage == "post_install" else 0,
                "",
                "setup failed" if failure_stage == "post_install" else "",
            )
        if command[-2:] == ["daemon", "restart"]:
            return plugin_manager.subprocess.CompletedProcess(
                command, 1, "", "restart failed"
            )
        return plugin_manager.subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(plugin_manager.subprocess, "run", run)
    plugin_config = {
        "title": "My Test Plugin",
        "description": "A test plugin",
        "pip": "my-plugin>=1.2.9",
    }
    if failure_stage == "post_install":
        plugin_config["post_install"] = "setup"
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {"my-plugin": plugin_config},
    )

    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]
    plugin._on_update(None)

    failure = plugin_state.get_activation_failure("my-plugin")
    assert failure == plugin.activation_error
    assert "daemon was not restarted" in failure or "daemon restart failed" in failure
    assert plugin.retry_activation_button.layout.display == ""
    assert plugin.version_warning.value
    assert (any(command[-2:] == ["daemon", "restart"] for command in commands)) == (
        failure_stage == "daemon_restart"
    )


def test_update_package_keeps_warning_if_minimum_is_not_met(monkeypatch):
    monkeypatch.setattr(
        QeAppPlugin,
        "_execute_command",
        lambda *_args, **_kwargs: True,
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.2.8", False, None),
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.2.8", False, None),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)
    subprocess_commands = []
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: (
            subprocess_commands.append(command)
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin>=1.2.9",
            }
        },
    )
    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]
    plugin._on_update(None)

    assert 'class="alert alert-danger"' in plugin.version_warning.value
    assert not plugin.update_button.disabled
    assert "does not meet the required version" in plugin.message_container.value
    assert subprocess_commands == [
        [
            plugin_manager.sys.executable,
            "-m",
            "aiidalab_qe",
            "test-plugin",
            "my-plugin",
        ]
    ]


def test_remove_reconciles_all_plugin_controls(monkeypatch):
    installed = {"version": "1.0"}
    commands = []
    daemon_commands = []

    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: (
            (installed["version"], False, None)
            if installed["version"]
            else (None, True, None)
        ),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)

    def execute(_row, command, *_args, **_kwargs):
        commands.append(command)
        installed["version"] = None
        return True

    monkeypatch.setattr(QeAppPlugin, "_execute_command", execute)
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: (
            daemon_commands.append(command)
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin>=2.0",
            }
        },
    )
    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]
    assert not plugin.update_button.disabled
    assert not plugin.remove_button.disabled

    plugin._on_remove(None)

    assert commands == [
        [plugin_manager.sys.executable, "-m", "pip", "uninstall", "-y", "my-plugin"]
    ]
    assert not plugin.install_button.disabled
    assert plugin.update_button.disabled
    assert plugin.remove_button.disabled
    assert plugin.version_warning.value == ""
    assert manager.accordion.get_title(0) == "My Test Plugin"
    assert daemon_commands == [["verdi", "daemon", "restart"]]


def test_install_reconciles_plugin_controls(monkeypatch):
    installed = {"version": "1.0"}
    commands = []
    daemon_commands = []

    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: (installed["version"], installed["version"] == "2.0", None),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)

    def execute(_row, command, *_args, **_kwargs):
        commands.append(command)
        installed["version"] = "2.0"
        return True

    monkeypatch.setattr(QeAppPlugin, "_execute_command", execute)

    class Result:
        returncode = 0
        stdout = ""
        stderr = ""

    def run(command, **_kwargs):
        daemon_commands.append(command)
        return Result()

    monkeypatch.setattr(plugin_manager.subprocess, "run", run)
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin>=2.0",
            }
        },
    )
    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]

    plugin._on_install(None)

    assert commands == [
        [
            plugin_manager.sys.executable,
            "-m",
            "pip",
            "install",
            "my-plugin>=2.0",
            "--user",
        ]
    ]
    assert plugin.install_button.disabled
    assert plugin.update_button.disabled
    assert not plugin.remove_button.disabled
    assert plugin.version_warning.value == ""
    assert manager.accordion.get_title(0).endswith("✅")
    assert daemon_commands == [
        [
            plugin_manager.sys.executable,
            "-m",
            "aiidalab_qe",
            "test-plugin",
            "my-plugin",
        ],
        ["verdi", "daemon", "restart"],
    ]


def test_run_post_install_runs_only_configured_command(monkeypatch):
    calls = []

    class Result:
        returncode = 0
        stdout = "setup complete"
        stderr = ""

    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: calls.append(command) or Result(),
    )
    output = ipw.HTML()
    message = ipw.HTML()
    data = QeAppPluginData.from_mapping(
        "my-plugin",
        {
            "title": "My Test Plugin",
            "description": "A test plugin",
            "pip": "my-plugin",
            "post_install": "setup",
        },
    )
    plugin = QeAppPlugin(data, plugin_name="my-plugin")
    plugin.output_container = output
    plugin.message_container = message
    plugin._run_post_install()

    assert calls == [[plugin_manager.sys.executable, "-m", "my_plugin", "setup"]]
    assert "background-color: #3B3B3B" in output.value
    assert "color: #FFFFFF" in output.value
    assert "setup complete" in output.value
    assert "Post-install completed" in message.value


@pytest.mark.parametrize(
    ("returncode", "expected_failure"),
    [
        (0, None),
        (1, "Post-install failed; plugin setup may be incomplete."),
    ],
)
def test_standalone_post_install_updates_activation_state(
    monkeypatch, returncode, expected_failure
):
    plugin_state.set_activation_failure("my-plugin", "Previous activation failure")
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    monkeypatch.setattr(plugin_manager, "is_version_compatible", lambda *_args: True)
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: plugin_manager.subprocess.CompletedProcess(
            command, returncode, "setup output", ""
        ),
    )
    data = QeAppPluginData.from_mapping(
        "my-plugin",
        {
            "title": "My Test Plugin",
            "description": "A test plugin",
            "pip": "my-plugin",
            "post_install": "setup",
        },
    )
    plugin = QeAppPlugin(data, plugin_name="my-plugin")

    plugin._on_post_install(None)

    assert plugin_state.get_activation_failure("my-plugin") == expected_failure
    assert plugin.activation_error == expected_failure
    assert (plugin.retry_activation_button.layout.display == "") == bool(
        expected_failure
    )


def test_post_install_only_button_is_enabled_when_package_is_not_installed(
    monkeypatch,
):
    monkeypatch.setattr(
        PluginManager,
        "_load_config",
        lambda _self: {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
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

    post_install_button = manager.accordion.children[0].children[2].children[1]
    assert not post_install_button.disabled
