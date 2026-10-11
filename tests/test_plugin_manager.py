from io import StringIO

import ipywidgets as ipw
import pytest

from aiidalab_qe.app.utils import plugin_manager
from aiidalab_qe.app.utils.plugin_manager import (
    PluginManager,
    QeAppPlugin,
    QeAppPluginData,
)
from aiidalab_qe.plugins import registry as plugin_registry

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


def _plugin_registry(monkeypatch, plugins):
    monkeypatch.setattr(plugin_registry, "load_plugin_config", lambda: plugins)
    return plugin_registry.PluginRegistry.load()


def _set_plugin_manager_registry(monkeypatch, plugins):
    registry = plugin_registry.PluginRegistry(
        plugins={
            name: QeAppPluginData.from_mapping(name, data)
            for name, data in plugins.items()
        }
    )
    monkeypatch.setattr(plugin_manager, "get_plugin_registry", lambda: registry)


def test_plugin_manager_uses_packaged_registry():
    """
    Test that PluginManager loads the packaged YAML and builds its UI properly.
    """
    manager = PluginManager()
    manager._build_ui()
    assert manager.registry is plugin_registry.get_plugin_registry()
    assert len(manager.accordion.children) == len(manager.plugins)
    for index, plugin in enumerate(manager.plugins.values()):
        assert plugin.title in manager.accordion.get_title(index)


def test_checked_in_plugin_registry_is_valid():
    registry = plugin_registry.load_plugin_config()

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


def test_plugin_registry_load_reports_type_errors(monkeypatch):
    registry = _plugin_registry(monkeypatch, {"my-plugin": []})

    assert not registry.plugins
    assert "plugin entry must be a mapping" in registry.config_error


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


def test_empty_registry_is_normalized_to_mapping(monkeypatch):
    registry = _plugin_registry(monkeypatch, {})
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_registry",
        lambda: registry,
    )
    manager = PluginManager()

    assert manager.data == {}
    assert manager.config_error is None


def test_invalid_registry_entry_is_reported_in_store(monkeypatch):
    registry = _plugin_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My plugin",
                "description": "Test",
                "pip": "not a requirement",
            }
        },
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_registry",
        lambda: registry,
    )
    displayed = []
    monkeypatch.setattr(plugin_manager, "display", displayed.append)

    manager = PluginManager()
    manager.display_ui()

    assert "invalid 'pip' requirement" in manager.config_error
    assert "invalid &#x27;pip&#x27; requirement" in displayed[0].value


def test_outdated_plugin_actions(monkeypatch):
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.0", False, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

    registry = _plugin_registry(
        monkeypatch,
        plugin_registry.yaml.safe_load(yaml_content),
    )
    monkeypatch.setattr(
        plugin_manager,
        "get_plugin_registry",
        lambda: registry,
    )
    manager = PluginManager()
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
    _set_plugin_manager_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin",
            }
        },
    )
    monkeypatch.setattr(
        plugin_registry, "get_plugin_version_info", lambda *_args: (None, True, None)
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

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
    _set_plugin_manager_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin>=1",
            }
        },
    )
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.0", plugin_compatible, None),
    )
    monkeypatch.setattr(
        plugin_registry, "is_version_compatible", lambda *_args: app_compatible
    )
    monkeypatch.setattr(
        plugin_registry,
        "poll_plugin_registration",
        lambda *_args: None,
    )

    manager = PluginManager()
    manager._build_ui()

    assert manager.accordion.get_title(0).endswith(expected_icon)


def test_plugin_manager_displays_persistent_ui_with_refresh_button(monkeypatch):
    _set_plugin_manager_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin",
            }
        },
    )
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (None, True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    monkeypatch.setattr(
        plugin_registry, "poll_plugin_registration", lambda *_args: None
    )
    displayed = []
    monkeypatch.setattr(plugin_manager, "display", displayed.append)

    manager = PluginManager()
    loading_states = []
    original_build_ui = PluginManager._build_ui

    def build_ui(self):
        loading_states.append(self.accordion.children[0].value)
        original_build_ui(self)

    monkeypatch.setattr(PluginManager, "_build_ui", build_ui)
    manager.display_ui()
    first_row = manager.accordion.children[0]
    first_row.message_container.value = "Keep this log"
    manager.display_ui()

    assert loading_states == ["Loading/validating plugins..."]
    assert displayed == [manager.ui, manager.ui]
    assert manager.ui.children == (
        manager.refresh_button,
        manager.daemon_warning,
        manager.accordion,
    )
    assert manager.accordion.children[0] is first_row
    assert first_row.message_container.value == "Keep this log"


def test_manager_refresh_updates_existing_rows_without_clearing_logs(monkeypatch):
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
        plugin_registry,
        "get_plugin_version_info",
        lambda name, *_args: version_info[name],
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    monkeypatch.setattr(
        plugin_registry, "poll_plugin_registration", lambda *_args: None
    )
    monkeypatch.setattr(plugin_manager, "poll_daemon_status", lambda: None)

    manager = PluginManager()
    manager._build_ui()
    first_row, second_row = manager.accordion.children
    first_row.message_container.value = "Previous update log"
    first_row.output_container.value = "Previous command output"
    version_info["first-plugin"] = ("1.0", False, None)
    version_info["second-plugin"] = ("2.0", True, None)

    manager.refresh_button.click()

    assert manager.accordion.children == (first_row, second_row)
    assert first_row.installed_version == "1.0"
    assert not first_row.plugin_compatible
    assert first_row.registration_error is None
    assert "Installed version 1.0" in first_row.version_warning.value
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
        plugin_registry.metadata, "version", lambda _package: installed_version
    )

    version, compatible, error = plugin_registry.get_plugin_version_info(
        "my-plugin", requirement
    )

    assert version == installed_version
    assert compatible is expected
    assert error is None


def test_get_plugin_version_info_not_installed(monkeypatch):
    def missing_package(_package):
        raise plugin_registry.metadata.PackageNotFoundError

    monkeypatch.setattr(plugin_registry.metadata, "version", missing_package)

    assert plugin_registry.get_plugin_version_info("my-plugin", "my-plugin>=1") == (
        None,
        True,
        None,
    )


def test_get_plugin_version_info_invalid_requirement():
    version, compatible, error = plugin_registry.get_plugin_version_info(
        "my-plugin", "not a requirement"
    )

    assert version is None
    assert not compatible
    assert error


def test_poll_post_install_skips_plugins_without_a_hook(monkeypatch):
    monkeypatch.setattr(
        plugin_registry.subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("no verifier should be run"),
    )

    assert plugin_registry.poll_post_install("my-plugin", None) == (None, None)


@pytest.mark.parametrize(
    ("result", "expected_error"),
    [
        ((True, ""), None),
        ((False, "broken entry point"), "broken entry point"),
        (
            (False, ""),
            "Plugin entry-point validation failed without an error message.",
        ),
    ],
)
def test_poll_plugin_registration_normalizes_result(
    monkeypatch, result, expected_error
):
    from aiidalab_qe.plugins import utils as plugin_utils

    monkeypatch.setattr(plugin_utils, "test_plugin_functionality", lambda _name: result)

    error = plugin_registry.poll_plugin_registration("my-plugin")

    assert error == expected_error


def test_poll_post_install_runs_plugin_verifier(monkeypatch):
    calls = []
    monkeypatch.setattr(
        plugin_registry.subprocess,
        "run",
        lambda command, **kwargs: (
            calls.append((command, kwargs))
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    assert plugin_registry.poll_post_install("aiida-bader", "post-install", "test") == (
        True,
        None,
    )
    assert calls[0][0] == [
        plugin_manager.sys.executable,
        "-m",
        "aiida_bader",
        "verify-post-install",
        "--profile",
        "test",
    ]
    assert calls[0][1]["timeout"] == plugin_registry.POST_INSTALL_VERIFY_TIMEOUT


def test_poll_post_install_reports_verifier_failure(monkeypatch):
    monkeypatch.setattr(
        plugin_registry.subprocess,
        "run",
        lambda command, **_kwargs: plugin_manager.subprocess.CompletedProcess(
            command, 1, "", "missing AiiDA code"
        ),
    )

    verified, error = plugin_registry.poll_post_install("my-plugin", "setup")

    assert verified is False
    assert "missing AiiDA code" in error


def test_poll_post_install_reports_timeout(monkeypatch):
    def timeout(*_args, **_kwargs):
        raise plugin_registry.subprocess.TimeoutExpired("verifier", 60)

    monkeypatch.setattr(plugin_registry.subprocess, "run", timeout)

    verified, error = plugin_registry.poll_post_install("my-plugin", "setup")

    assert verified is False
    assert "timed out" in error


def test_plugin_status_marks_failed_post_install_incompatible(monkeypatch):
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.0", True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    monkeypatch.setattr(
        plugin_registry,
        "poll_post_install",
        lambda *_args: (False, "Post-install check failed"),
    )
    plugin_data = QeAppPluginData.from_mapping(
        "my-plugin",
        {
            "title": "My plugin",
            "description": "Test plugin",
            "pip": "my-plugin",
            "post_install": "setup",
        },
    )

    status = plugin_registry.get_plugin_status("my-plugin", plugin_data)

    assert status.post_install_verified is False
    assert status.post_install_error == "Post-install check failed"
    assert status.is_incompatible


def test_update_package_clears_warning_after_minimum_is_met(monkeypatch):
    commands = []

    def execute(_row, command, *_args, **_kwargs):
        commands.append(command)
        installed_version["value"] = "1.2.9"
        return True

    installed_version = {"value": "1.2.8"}
    monkeypatch.setattr(QeAppPlugin, "_execute_command", execute)
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (
            installed_version["value"],
            installed_version["value"] == "1.2.9",
            None,
        ),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    daemon_commands = []
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: (
            daemon_commands.append(command)
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
        plugin_registry.metadata,
        "version",
        lambda _name: installed_version["value"],
    )
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (
            installed_version["value"],
            installed_version["value"] == "1.2.9",
            None,
        ),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    subprocess_commands = []
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: (
            subprocess_commands.append(command)
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (
            installed_version["value"],
            installed_version["value"] == "1.2.9",
            None,
        ),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

    subprocess_commands = []

    def fail_plugin_test(command, **_kwargs):
        subprocess_commands.append(command)
        return plugin_manager.subprocess.CompletedProcess(command, 1, "", "failed")

    monkeypatch.setattr(plugin_manager.subprocess, "run", fail_plugin_test)
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
    assert plugin.registration_error
    assert plugin.retry_activation_button.layout.display == ""
    assert plugin.update_button.disabled is False


def test_retry_activation_clears_failure_and_restarts_daemon(monkeypatch):
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    monkeypatch.setattr(
        plugin_registry,
        "poll_plugin_registration",
        lambda *_args: "Previous validation failed",
    )
    calls = []

    def run(command, **_kwargs):
        calls.append(command)
        return plugin_manager.subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(plugin_manager.subprocess, "run", run)
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
    assert plugin.registration_error is None
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
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
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
    _set_plugin_manager_registry(monkeypatch, {"my-plugin": plugin_config})

    manager = PluginManager()
    manager._build_ui()
    plugin = manager.accordion.children[0]
    plugin._on_update(None)

    if failure_stage == "post_install":
        assert plugin.post_install_error
    else:
        assert plugin.daemon_error
    assert plugin.retry_activation_button.layout.display == ""
    assert (
        plugin.version_warning.value
        if failure_stage == "post_install"
        else plugin.daemon_error
    )
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
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.2.8", False, None),
    )
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.2.8", False, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    subprocess_commands = []
    monkeypatch.setattr(
        plugin_manager.subprocess,
        "run",
        lambda command, **_kwargs: (
            subprocess_commands.append(command)
            or plugin_manager.subprocess.CompletedProcess(command, 0, "", "")
        ),
    )
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (
            (installed["version"], False, None)
            if installed["version"]
            else (None, True, None)
        ),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

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
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: (installed["version"], installed["version"] == "2.0", None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

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
    _set_plugin_manager_registry(
        monkeypatch,
        {
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
    ("returncode", "expected_error"),
    [
        (0, None),
        (1, "Post-install verifier failed"),
    ],
)
def test_standalone_post_install_refreshes_runtime_status(
    monkeypatch, returncode, expected_error
):
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    monkeypatch.setattr(
        plugin_registry,
        "poll_plugin_registration",
        lambda *_args: None,
    )
    monkeypatch.setattr(
        plugin_registry,
        "poll_post_install",
        lambda *_args: (returncode == 0, expected_error),
    )
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

    assert plugin.post_install_error == expected_error
    assert (plugin.retry_activation_button.layout.display == "") == bool(expected_error)


def test_runtime_post_install_failure_warns_and_offers_retry(monkeypatch):
    monkeypatch.setattr(
        plugin_registry,
        "get_plugin_version_info",
        lambda *_args: ("1.2.9", True, None),
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)
    monkeypatch.setattr(
        plugin_registry,
        "poll_post_install",
        lambda *_args: (False, "AiiDA code phonopy@localhost is missing."),
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

    assert "phonopy@localhost is missing" in plugin.version_warning.value
    assert plugin.retry_activation_button.layout.display == ""


def test_post_install_only_button_is_enabled_when_package_is_not_installed(
    monkeypatch,
):
    _set_plugin_manager_registry(
        monkeypatch,
        {
            "my-plugin": {
                "title": "My Test Plugin",
                "description": "A test plugin",
                "pip": "my-plugin",
                "post_install": "setup",
            }
        },
    )
    monkeypatch.setattr(
        plugin_registry, "get_plugin_version_info", lambda *_args: (None, True, None)
    )
    monkeypatch.setattr(plugin_registry, "is_version_compatible", lambda *_args: True)

    manager = PluginManager()
    manager._build_ui()

    post_install_button = manager.accordion.children[0].children[2].children[1]
    assert not post_install_button.disabled
