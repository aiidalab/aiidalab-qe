"""
plugin_manager.py

Module that provides a plugin management interface, allowing the user to
install and remove Python-based AiiDAlab plugins, with real-time streaming
of command output in a Jupyter environment.
"""

from __future__ import annotations

import html
import logging
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

import ipywidgets as ipw
import traitlets as tl
from IPython.display import display
from packaging.requirements import Requirement

from aiidalab_qe import __version__
from aiidalab_qe.plugins.registry import (
    DEFAULT_PLUGIN_CONFIG_SOURCE,
    PluginRegistry,
    QeAppPluginData,
    get_plugin_version_info,
    is_version_compatible,
    load_plugin_config,
)
from aiidalab_qe.plugins.state import (
    clear_activation_failure,
    get_activation_failure,
    set_activation_failure,
)

LOGGER = logging.getLogger(__name__)

COLOR_MAP = {
    "experimental": "#FF8C00",  # 🟠 Orange - Early development
    "beta": "#FFA500",  # 🟡 Darker Orange - Testing phase
    "stable": "#1E90FF",  # 🔵 Dodger Blue - Well-tested, feature-complete
    "production": "#008000",  # 🟢 Green - Fully stable and recommended
    "deprecated": "#FF0000",  # 🔴 Red - No longer maintained
    "archived": "#808080",  # ⚪ Grey - Retained for reference, no updates
}

BUTTON_WIDTH = "120px"


class QeAppPlugin(ipw.VBox):
    """Widget and actions for one validated plugin registry entry."""

    installed_version = tl.Unicode(default_value=None, allow_none=True)
    plugin_compatible = tl.Bool(default_value=True)
    app_compatible = tl.Bool(default_value=True)
    requirement_error = tl.Unicode(default_value=None, allow_none=True)
    activation_error = tl.Unicode(default_value=None, allow_none=True)
    version_warning_text = tl.Unicode(default_value="")

    def __init__(
        self,
        data: QeAppPluginData,
        plugin_name: str,
        title_updater: Callable[[str], None] | None = None,
        **kwargs,
    ):
        self.data = data
        self.plugin_name = plugin_name
        self.title_updater = title_updater

        self.version_warning = ipw.HTML()
        tl.dlink(
            (self, "version_warning_text"),
            (self.version_warning, "value"),
        )

        self.message_container = ipw.HTML(
            value="",
            layout=ipw.Layout(
                max_height="250px",
                overflow="auto",
                border="1px solid #9e9e9e",
                display="none",
            ),
        )
        self.output_container = ipw.HTML(
            value="",
            layout=ipw.Layout(
                max_height="250px",
                overflow="auto",
                border="1px solid #9e9e9e",
                display="none",
            ),
        )
        self.install_button = ipw.Button(
            description="Install",
            tooltip="Install the plugin",
            button_style="success",
            layout=ipw.Layout(width=BUTTON_WIDTH),
        )
        self.install_button.on_click(self._on_install)

        self.post_install_button = ipw.Button(
            description="Run post-install",
            tooltip="Run the post-install script",
            button_style="info",
            layout=ipw.Layout(
                width=BUTTON_WIDTH,
                display="" if self.data.post_install else "none",
            ),
        )
        self.post_install_button.on_click(self._on_post_install)

        self.update_button = ipw.Button(
            description="Update",
            tooltip="Update the plugin",
            button_style="warning",
            layout=ipw.Layout(width=BUTTON_WIDTH),
        )
        self.update_button.on_click(self._on_update)

        self.remove_button = ipw.Button(
            description="Remove",
            tooltip="Remove the plugin",
            button_style="danger",
            layout=ipw.Layout(width=BUTTON_WIDTH),
        )
        self.remove_button.on_click(self._on_remove)

        self.retry_activation_button = ipw.Button(
            description="Retry",
            tooltip="Retry plugin activation",
            button_style="info",
            layout=ipw.Layout(width=BUTTON_WIDTH, display="none"),
        )
        self.retry_activation_button.on_click(self._on_retry_activation)

        self.clear_output_button = ipw.Button(
            description="Clear output",
            tooltip="Clear plugin messages and command output",
            disabled=True,
            layout=ipw.Layout(width=BUTTON_WIDTH),
        )
        self.clear_output_button.on_click(self._on_clear_output)

        for output_widget in (self.message_container, self.output_container):
            output_widget.observe(
                self._sync_clear_output_button,
                "value",
            )

        super().__init__(
            children=[
                ipw.HTML(self._details_html()),
                self.version_warning,
                ipw.HBox(
                    children=[
                        self.install_button,
                        self.post_install_button,
                        self.update_button,
                        self.remove_button,
                        self.retry_activation_button,
                        self.clear_output_button,
                    ]
                ),
                self.message_container,
                self.output_container,
            ],
            **kwargs,
        )

        # This separate call ensures that the UI is refreshed even if
        # the state remains the same after the plugin state refresh.
        self._refresh_ui_state()

        self.refresh_plugin_state()

    @property
    def is_installed(self) -> bool:
        return self.installed_version is not None

    def refresh_plugin_state(self) -> None:
        """Refresh plugin state from package metadata and persisted activation status."""
        installed_version, plugin_compatible, requirement_error = (
            get_plugin_version_info(self.plugin_name, self.data.pip or self.plugin_name)
        )
        app_compatible = is_version_compatible(self.data.requires_aiidalab_qe)
        activation_error = get_activation_failure(self.plugin_name)

        with self.hold_trait_notifications():
            self.installed_version = installed_version
            self.plugin_compatible = plugin_compatible
            self.requirement_error = requirement_error
            self.app_compatible = app_compatible
            self.activation_error = activation_error

    def update_title(self) -> None:
        if self.title_updater is None:
            return

        status = "⚠️" if self.version_warning_text else "✅" if self.is_installed else ""

        self.title_updater(f"{self.data.title} {status}".rstrip())

    @tl.observe(
        "installed_version",
        "plugin_compatible",
        "app_compatible",
        "requirement_error",
        "activation_error",
    )
    def _refresh_ui_state(self, _change: dict | None = None) -> None:
        """Update the UI state of the plugin."""
        self._update_version_warning()
        self._update_action_buttons()
        self.update_title()

    def _update_version_warning(self) -> None:
        """Build and store the warning text from the current plugin state."""
        warning = ""

        if self.activation_error:
            warning += (
                '<div class="alert alert-danger" role="alert">'
                "Plugin activation is not confirmed: "
                f"{html.escape(self.activation_error)}</div>"
            )

        if not self.app_compatible:
            warning += (
                '<div class="alert alert-danger" role="alert">'
                f"This plugin requires aiidalab_qe {self.data.requires_aiidalab_qe}, "
                f"but you have {__version__}.</div>"
            )

        if self.requirement_error:
            warning += (
                '<div class="alert alert-danger" role="alert">'
                "Could not evaluate plugin requirement/version: "
                f"{self.requirement_error}</div>"
            )
        elif self.is_installed and not self.plugin_compatible:
            requirement = Requirement(self.data.pip or self.plugin_name)
            warning += (
                '<div class="alert alert-danger" role="alert">'
                f"Installed version {self.installed_version} does not satisfy "
                f"{requirement.name}{requirement.specifier}. Please update accordingly."
                "</div>"
            )

        self.version_warning_text = warning

    def _update_action_buttons(self) -> None:
        """Set action availability from the current plugin state."""
        if self.install_button is not None:
            self.install_button.disabled = (
                self.is_installed
                or not self.app_compatible
                or bool(self.requirement_error)
            )

        if self.update_button is not None:
            self.update_button.disabled = (
                not self.is_installed
                or (self.plugin_compatible and not self.activation_error)
                or not self.app_compatible
            )

        if self.remove_button is not None:
            self.remove_button.disabled = not self.is_installed

        if self.retry_activation_button is not None:
            self.retry_activation_button.layout.display = (
                "" if self.activation_error and self.is_installed else "none"
            )
            self.retry_activation_button.disabled = not self.app_compatible

    def _details_html(self) -> str:
        status_value = self.data.status.strip().lower()
        badge_color = COLOR_MAP.get(status_value, "#666666")
        display_text = status_value.capitalize() if status_value else "N/A"

        badge_html = (
            f'<span style="background-color: {badge_color}; color: #FFFFFF; '
            f'border-radius: 4px; padding: 2px 6px;">{display_text}</span>'
        )

        details = (
            f"<b>Package:</b> {self.plugin_name}<br>"
            f"<b>Author:</b> {self.data.author}<br>"
            f"<b>Description:</b> {self.data.description}<br>"
            f"<b>Status:</b> {badge_html}<br>"
        )

        if self.data.documentation:
            details += (
                f"<b>Documentation:</b> <a href='{self.data.documentation}' "
                "target='_blank'>Visit</a><br>"
            )

        if self.data.github:
            details += (
                f"<b>Github:</b> <a href='{self.data.github}' target='_blank'>Visit</a>"
            )

        return details

    def _on_install(self, _button: ipw.Button) -> None:
        message = self.message_container
        output = self.output_container
        output.layout.display = "block"
        message.layout.display = "block"
        message.value = ""
        output.value = ""

        self._append_message(f"Installing {self.plugin_name}...")
        set_activation_failure(
            self.plugin_name,
            "Installation is in progress; activation has not been confirmed.",
        )

        install_source = self.data.pip or f"git+{self.data.github}"
        installed = self._execute_command(
            [sys.executable, "-m", "pip", "install", install_source, "--user"],
        )
        if not installed:
            set_activation_failure(
                self.plugin_name,
                "Installation did not complete; the plugin state may have changed.",
            )
            self._append_message(
                "Installation did not complete. The plugin state may have changed, "
                "so activation is not confirmed. Review the command output.",
                color="#FF0000",
            )
            self.refresh_plugin_state()
            return

        if not self._validate_plugin_installation(remove_on_test_failure=True):
            self.refresh_plugin_state()
            if self.is_installed:
                set_activation_failure(
                    self.plugin_name,
                    "Installation setup or plugin validation failed. The plugin "
                    "remains installed, and the daemon was not restarted.",
                )
            else:
                clear_activation_failure(self.plugin_name)
            self.refresh_plugin_state()
            return

        if not self._restart_daemon():
            set_activation_failure(
                self.plugin_name,
                "Plugin validation passed, but the daemon restart failed.",
            )
            self._append_message(
                "Plugin validation passed, but the daemon restart failed. "
                "Activation is not confirmed.",
                color="#FF0000",
            )
            self.refresh_plugin_state()
            return

        clear_activation_failure(self.plugin_name)
        self._append_message("Plugin installed successfully.", color="#008000")

        self.refresh_plugin_state()

    def _validate_plugin_installation(
        self, remove_on_test_failure: bool = False
    ) -> bool:
        """Run the plugin setup hook and verify that its entry points load."""
        if self.data.post_install and not self._run_post_install(clear_output=False):
            return False

        self._append_message("Testing plugin loading...", color="#008000")

        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "aiidalab_qe",
                    "test-plugin",
                    self.plugin_name,
                ],
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError:
            LOGGER.exception(
                "Could not run plugin loading test for %s", self.plugin_name
            )
            self._append_message(
                "The plugin test could not be started. Details were logged for debugging.",
                color="#FF0000",
            )
            return False

        if result.stdout:
            self._append_output(result.stdout)
        if result.stderr:
            self._append_output(result.stderr)

        if result.returncode == 0:
            self._append_message("Plugin test passed.", color="#008000")
            return True

        LOGGER.error(
            "Plugin test failed for %s (exit code %s): %s%s",
            self.plugin_name,
            result.returncode,
            result.stdout,
            result.stderr,
        )
        message = f"The plugin test for {self.plugin_name} did not pass."
        if remove_on_test_failure:
            message += (
                " The plugin will be removed to prevent use in an incomplete state."
            )
        else:
            message += (
                " The plugin remains installed, but activation is not confirmed. "
                "The daemon was not restarted."
            )
        self._append_message(message + " Details were logged.", color="#FF0000")

        if remove_on_test_failure:
            self._remove_package(clear_output=False, restart_daemon=False)

        return False

    def _on_update(self, _button: ipw.Button) -> None:
        output = self.output_container
        message = self.message_container
        output.layout.display = "block"
        message.layout.display = "block"
        message.value = ""

        self._append_message(f"Updating {self.plugin_name}...")
        set_activation_failure(
            self.plugin_name,
            "Update is in progress; activation has not been confirmed.",
        )

        requirement = self.data.pip or f"git+{self.data.github}"
        result = self._execute_command(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--upgrade",
                requirement,
                "--user",
            ]
        )
        if not result:
            set_activation_failure(
                self.plugin_name,
                "Update did not complete; the plugin state may have changed.",
            )
            self._append_message(
                "Update did not complete. The plugin state may have changed, "
                "so activation is not confirmed and the daemon was not restarted.",
                color="#FF0000",
            )
            self.refresh_plugin_state()
            return

        if not self._validate_plugin_installation():
            set_activation_failure(
                self.plugin_name,
                "Setup or plugin validation failed after the plugin update. "
                "The plugin remains installed, but the daemon was not restarted.",
            )
            self.refresh_plugin_state()
            return

        self.refresh_plugin_state()

        if (
            self.is_installed
            and self.plugin_compatible
            and self.app_compatible
            and not self.requirement_error
        ):
            if not self._restart_daemon():
                set_activation_failure(
                    self.plugin_name,
                    "Plugin validation passed, but the daemon restart failed.",
                )
                self._append_message(
                    "Plugin validation passed, but the daemon restart failed. "
                    "Activation is not confirmed.",
                    color="#FF0000",
                )
                self.refresh_plugin_state()
                return
            clear_activation_failure(self.plugin_name)
            self.refresh_plugin_state()
            self._append_message(
                f"Updated {self.plugin_name} to {self.installed_version}.",
                color="#008000",
            )
        elif result:
            reason = self.requirement_error or (
                f"Installed version {self.installed_version or 'unknown'} does not meet "
                "the required version."
            )
            set_activation_failure(
                self.plugin_name,
                f"Update completed, but activation was not confirmed: {reason}",
            )
            self.refresh_plugin_state()
            self._append_message(reason, color="#FF0000")

    def _on_retry_activation(self, _button: ipw.Button) -> None:
        self.message_container.layout.display = "block"
        self.output_container.layout.display = "block"
        self.message_container.value = ""
        set_activation_failure(
            self.plugin_name,
            "Activation retry is in progress; activation has not been confirmed.",
        )

        if not self._validate_plugin_installation():
            set_activation_failure(
                self.plugin_name,
                "Setup or plugin validation failed. The daemon was not restarted.",
            )
            self.refresh_plugin_state()
            return

        self.refresh_plugin_state()
        if (
            not self.plugin_compatible
            or not self.app_compatible
            or self.requirement_error
        ):
            set_activation_failure(
                self.plugin_name,
                "Validation passed, but the installed plugin does not satisfy the "
                "registered compatibility requirements.",
            )
            self.refresh_plugin_state()
            return

        if not self._restart_daemon():
            set_activation_failure(
                self.plugin_name,
                "Plugin validation passed, but the daemon restart failed.",
            )
            self._append_message(
                "Plugin validation passed, but the daemon restart failed. "
                "Activation is not confirmed.",
                color="#FF0000",
            )
            self.refresh_plugin_state()
            return

        clear_activation_failure(self.plugin_name)
        self._append_message("Plugin activation succeeded.", color="#008000")
        self.refresh_plugin_state()

    def _on_remove(self, _button: ipw.Button) -> None:
        self.message_container.layout.display = "block"
        self.output_container.layout.display = "block"
        self._remove_package()

    def _remove_package(
        self, clear_output: bool = True, restart_daemon: bool = True
    ) -> bool:
        self._append_message(f"Removing {self.plugin_name}...")

        requirement = Requirement(self.data.pip or self.plugin_name)
        result = self._execute_command(
            [sys.executable, "-m", "pip", "uninstall", "-y", requirement.name],
            clear_output=clear_output,
        )
        if result:
            restart_succeeded = not restart_daemon or self._restart_daemon()
            if restart_succeeded:
                clear_activation_failure(self.plugin_name)
            else:
                set_activation_failure(
                    self.plugin_name,
                    "Package removal succeeded, but the daemon restart failed.",
                )
            self._append_message(
                f"{self.plugin_name} removed successfully.",
                color="#008000",
            )
            if not restart_succeeded:
                self._append_message(
                    "The daemon restart failed; reload the app before using plugins.",
                    color="#FF0000",
                )

        self.refresh_plugin_state()

        return result

    def _on_post_install(self, _button: ipw.Button) -> None:
        if self._run_post_install():
            clear_activation_failure(self.plugin_name)
        else:
            set_activation_failure(
                self.plugin_name,
                "Post-install failed; plugin setup may be incomplete.",
            )
        self.refresh_plugin_state()

    def _run_post_install(self, clear_output: bool = True) -> bool:
        output = self.output_container
        message = self.message_container
        output.layout.display = "block"
        message.layout.display = "block"

        if clear_output:
            output.value = ""
            message.value = ""

        self._append_message(f"Running post-install for {self.plugin_name}...")

        command = [
            sys.executable,
            "-m",
            self.plugin_name.replace("-", "_"),
            self.data.post_install,
        ]

        try:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        except OSError:
            LOGGER.exception("Could not start post-install for %s", self.plugin_name)
            self._append_message(
                "Post-install could not be started. Details were logged for debugging.",
                color="#FF0000",
            )
            return False

        if result.stdout:
            self._append_output(result.stdout)

        if result.returncode == 0:
            self._append_message(
                f"Post-install completed for {self.plugin_name}.",
                color="#008000",
            )
            return True

        LOGGER.error(
            "Post-install for %s failed with exit code %s: %s",
            self.plugin_name,
            result.returncode,
            result.stdout,
        )

        self._append_message(
            "Post-install did not complete. The plugin is present, but setup may be "
            "incomplete. Review the command output; details were logged for debugging.",
            color="#FF0000",
        )

        return False

    def _on_clear_output(self, _button: ipw.Button) -> None:
        self._clear_output()

    def _append_output(self, text: str) -> None:
        """Append escaped command output while preserving terminal line breaks."""
        self.output_container.value += (
            '<div style="background-color: #3B3B3B; color: #FFFFFF; '
            'white-space: pre-wrap; margin: 0; padding: 0 4px;">'
            + html.escape(text)
            + "</div>"
        )

    def _append_message(self, text: str, color: str = "#000000") -> None:
        """Append an escaped high-level status message to this row."""
        self.message_container.value += (
            f'<div style="color: {color}; padding: 0 4px;">{html.escape(text)}</div>'
        )

    def _execute_command(self, command: list, clear_output: bool = True) -> bool:
        """Execute a command and stream its output into this row."""
        if clear_output:
            self.output_container.value = ""
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except OSError:
            LOGGER.exception("Could not start plugin command %r", command)
            self._append_output(
                "The command could not be started. Details were logged for debugging.\n"
            )
            return False

        for output in process.stdout:
            self._append_output(output)

        return_code = process.wait()

        if return_code != 0:
            LOGGER.error(
                "Plugin command %r failed with exit code %s", command, return_code
            )
            self._append_output(
                "\nThe command did not complete successfully. Review the output above; "
                "details were logged for debugging.\n"
            )
            return False

        return True

    def _clear_output(self) -> None:
        """Clear and hide this row's messages and command output."""
        self.message_container.value = ""
        self.output_container.value = ""
        self.message_container.layout.display = "none"
        self.output_container.layout.display = "none"

    def _sync_clear_output_button(self, _change: dict) -> None:
        self.clear_output_button.disabled = not (
            self.message_container.value or self.output_container.value
        )

    @staticmethod
    def _restart_daemon() -> bool:
        try:
            result = subprocess.run(
                ["verdi", "daemon", "restart"], capture_output=True, check=False
            )
        except OSError:
            LOGGER.exception("Could not restart the AiiDA daemon")
            return False
        if result.returncode:
            LOGGER.error("Could not restart the AiiDA daemon: %s", result.stderr)
            return False
        return True


class PluginManager:
    """Plugin manager for AiiDAlab Quantum ESPRESSO.

    A manager class that reads a plugin configuration file (YAML),
    and creates an interactive Accordion UI for installing/uninstalling
    those plugins in a Jupyter environment.
    """

    def __init__(self, config_source: str | Path = DEFAULT_PLUGIN_CONFIG_SOURCE):
        """Initialize the PluginManager with a path to a YAML config or a URL.

        :param config_source: A local YAML path, a plugin resource, or a remote URL.
        """
        self.config_source = config_source
        self.registry = PluginRegistry.load(config_source, loader=self._load_config)
        self.config_error = self.registry.config_error
        self.plugins = self.registry.plugins

        self.accordion = ipw.Accordion()

        self.refresh_button = ipw.Button(
            description="Refresh",
            tooltip="Refresh plugin status",
            icon="refresh",
            button_style="primary",
            layout=ipw.Layout(width="fit-content", margin="0 2px 10px"),
        )
        self.refresh_button.on_click(self._on_refresh)

        self.ui = ipw.VBox(
            children=[
                self.refresh_button,
                self.accordion,
            ],
        )

        self._accordion_built = False

    @property
    def data(self) -> dict:
        """Return normalized registry entries for existing consumers."""
        return self.registry.data

    def display_ui(self) -> None:
        """Display the plugin management UI.

        Renders a warning if the plugin registry could not be loaded.
        """
        if self.config_error:
            display(
                ipw.HTML(
                    value=(
                        '<div class="alert alert-danger" role="alert">'
                        "Plugin registry could not be loaded: "
                        f"{html.escape(self.config_error)}</div>"
                    )
                )
            )
            return

        self._build_ui()

        display(self.ui)

    def refresh_plugins(self) -> None:
        """Refresh external status for each existing plugin row."""
        for plugin in self.accordion.children:
            plugin.refresh_plugin_state()

    def _on_refresh(self, _button: ipw.Button) -> None:
        self.refresh_plugins()

    def _load_config(self) -> dict:
        """Load raw YAML data from the configured source."""
        return load_plugin_config(self.config_source)

    def _build_ui(self) -> None:
        """Build the Accordion UI from the validated plugin records.

        Title updates are delegated to each individual plugin row.
        """
        if self._accordion_built:
            return

        def update_title(index: int, title: str) -> None:
            if index < len(self.accordion.children):
                self.accordion.set_title(index, title)

        self.accordion.children = [
            QeAppPlugin(
                plugin_data,
                plugin_name,
                title_updater=lambda title, index=index: update_title(index, title),
            )
            for index, (plugin_name, plugin_data) in enumerate(self.plugins.items())
        ]

        for plugin in self.accordion.children:
            plugin.update_title()

        self._accordion_built = True
