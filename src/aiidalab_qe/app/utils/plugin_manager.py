"""
plugin_manager.py

Module that provides a plugin management interface, allowing the user to
install and remove Python-based AiiDAlab plugins, with real-time streaming
of command output in a Jupyter environment.
"""

import html
import logging
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib import metadata
from pathlib import Path
from urllib.parse import urlparse

import ipywidgets as ipw
import requests
import yaml
from IPython.display import display
from packaging.requirements import InvalidRequirement, Requirement
from packaging.specifiers import InvalidSpecifier, SpecifierSet
from packaging.version import InvalidVersion, Version

LOGGER = logging.getLogger(__name__)

# Define badge colors based on status
COLOR_MAP = {
    "experimental": "#FF8C00",  # 🟠 Orange - Early development
    "beta": "#FFA500",  # 🟡 Darker Orange - Testing phase
    "stable": "#1E90FF",  # 🔵 Dodger Blue - Well-tested, feature-complete
    "production": "#008000",  # 🟢 Green - Fully stable and recommended
    "deprecated": "#FF0000",  # 🔴 Red - No longer maintained
    "archived": "#808080",  # ⚪ Grey - Retained for reference, no updates
}
DEFAULT_PLUGIN_CONFIG_SOURCE = (
    "https://raw.githubusercontent.com/aiidalab/aiidalab-qe/main/plugins.yaml"
)
BUTTON_WIDTH = "120px"


def get_aiidalab_qe_version() -> str:
    """
    Get the installed version of aiidalab_qe.
    Returns 'unknown' if the package is not installed.
    """
    import aiidalab_qe

    try:
        return aiidalab_qe.__version__
    except AttributeError:
        return "unknown"


INSTALLED_AIIDA_QE_VERSION = get_aiidalab_qe_version()


def is_version_compatible(required_version: str) -> bool:
    """
    Check if the installed aiidalab_qe version satisfies the required version constraint.
    This function explicitly allows pre-release versions if they match the specifier.
    """
    if not required_version:
        return True

    if INSTALLED_AIIDA_QE_VERSION == "unknown":
        return False  # aiidalab_qe is not installed

    try:
        specifier = SpecifierSet(required_version)  # Handles constraints like '>=24.10'
        installed_version = Version(INSTALLED_AIIDA_QE_VERSION)
        return specifier.contains(installed_version, prereleases=True)

    except (InvalidSpecifier, InvalidVersion) as error:
        LOGGER.warning(
            "Could not parse version requirement %r: %s",
            required_version,
            error,
        )
        return False


def is_package_installed(package_name: str) -> bool:
    """
    Check if a given Python package is already installed.
    """
    try:
        metadata.version(package_name)
    except metadata.PackageNotFoundError:
        return False
    else:
        return True


def get_plugin_version_info(
    plugin_name: str,
    pip_requirement: str | None,
) -> tuple[str | None, bool, str | None]:
    """Return the installed version, compatibility, and any requirement error."""
    try:
        requirement = Requirement(pip_requirement or plugin_name)
    except InvalidRequirement as error:
        return None, False, str(error)

    try:
        installed_version = metadata.version(requirement.name)
    except metadata.PackageNotFoundError:
        return None, True, None

    try:
        version = Version(installed_version)
    except InvalidVersion as error:
        return installed_version, False, str(error)

    compatible = not requirement.specifier or requirement.specifier.contains(
        version,
        prereleases=True,
    )
    return installed_version, compatible, None


@dataclass(frozen=True)
class QeAppPluginData:
    registry_name: str
    title: str
    description: str
    package: str
    pip: str | None = None
    github: str | None = None
    author: str = "N/A"
    documentation: str | None = None
    post_install: str | None = None
    status: str = ""
    category: str = "calculation"
    requires_aiidalab_qe: str | None = None
    extra: dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_mapping(cls, registry_name: str, data: Mapping) -> "QeAppPluginData":
        """Create a plugin record after validating fields consumed by the app."""
        if not isinstance(registry_name, str) or not registry_name.strip():
            raise ValueError("plugin registry key must be a non-empty string")
        if not isinstance(data, Mapping):
            raise TypeError(f"{registry_name}: plugin entry must be a mapping")

        def required_string(field: str) -> str:
            value = data.get(field)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"{registry_name}: '{field}' must be a non-empty string"
                )
            return value

        def optional_string(field: str, default: str | None = None) -> str | None:
            value = data.get(field, default)
            if value is None:
                return None
            if not isinstance(value, str):
                raise TypeError(f"{registry_name}: '{field}' must be a string")
            return value or None

        title = required_string("title")
        description = required_string("description")

        pip_requirement = optional_string("pip")
        github = optional_string("github")
        if not pip_requirement and not github:
            raise ValueError(f"{registry_name}: either 'pip' or 'github' is required")

        package = optional_string("package", registry_name) or registry_name

        try:
            Requirement(pip_requirement or package)
        except InvalidRequirement as error:
            field_name = "pip" if pip_requirement else "package"
            raise ValueError(
                f"{registry_name}: invalid '{field_name}' requirement: {error}"
            ) from error

        app_requirement = optional_string("requires_aiidalab_qe")
        if app_requirement:
            try:
                SpecifierSet(app_requirement)
            except InvalidSpecifier as error:
                raise ValueError(
                    f"{registry_name}: invalid 'requires_aiidalab_qe' specifier: {error}"
                ) from error

        known_fields = {
            "title",
            "description",
            "package",
            "pip",
            "github",
            "author",
            "documentation",
            "post_install",
            "status",
            "category",
            "requires_aiidalab_qe",
        }

        return cls(
            registry_name=registry_name,
            title=title,
            description=description,
            package=package,
            pip=pip_requirement,
            github=github,
            author=optional_string("author", "N/A") or "N/A",
            documentation=optional_string("documentation"),
            post_install=optional_string("post_install"),
            status=optional_string("status", "") or "",
            category=optional_string("category", "calculation") or "calculation",
            requires_aiidalab_qe=app_requirement,
            extra={
                key: value for key, value in data.items() if key not in known_fields
            },
        )

    def as_mapping(self) -> dict:
        """Expose normalized metadata for existing registry consumers."""
        return {
            **self.extra,
            "title": self.title,
            "description": self.description,
            "package": self.package,
            "pip": self.pip,
            "github": self.github,
            "author": self.author,
            "documentation": self.documentation,
            "post_install": self.post_install,
            "status": self.status,
            "category": self.category,
            "requires_aiidalab_qe": self.requires_aiidalab_qe,
        }


class QeAppPluginRow:
    """Widgets and actions for one validated plugin registry entry."""

    def __init__(self, data: QeAppPluginData):
        self.data = data
        self.accordion: ipw.Accordion | None = None
        self.index: int | None = None
        self.installed_version: str | None = None
        self.plugin_compatible = True
        self.app_compatible = True
        self.requirement_error: str | None = None
        self.version_warning: ipw.HTML | None = None
        self.message_container: ipw.HTML | None = None
        self.output_container: ipw.HTML | None = None
        self.install_button: ipw.Button | None = None
        self.post_install_button: ipw.Button | None = None
        self.update_button: ipw.Button | None = None
        self.remove_button: ipw.Button | None = None
        self.clear_output_button: ipw.Button | None = None

    @property
    def is_installed(self) -> bool:
        return self.installed_version is not None

    def build_panel(self, accordion: ipw.Accordion, index: int) -> ipw.VBox:
        self.accordion = accordion
        self.index = index
        self.version_warning = ipw.HTML()
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
            button_style="success",
            layout=ipw.Layout(width=BUTTON_WIDTH),
        )
        self.install_button.on_click(self._on_install)

        self.post_install_button = ipw.Button(
            description="Run post-install",
            button_style="info",
            layout=ipw.Layout(
                width=BUTTON_WIDTH,
                display="" if self.data.post_install else "none",
            ),
        )
        self.post_install_button.on_click(self._on_post_install)

        self.update_button = ipw.Button(
            description="Update",
            button_style="warning",
            layout=ipw.Layout(width=BUTTON_WIDTH),
        )
        self.update_button.on_click(self._on_update)

        self.remove_button = ipw.Button(
            description="Remove",
            button_style="danger",
            layout=ipw.Layout(width=BUTTON_WIDTH),
        )
        self.remove_button.on_click(self._on_remove)

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

        self.reconcile_state()

        return ipw.VBox(
            children=[
                ipw.HTML(self._details_html()),
                self.version_warning,
                ipw.HBox(
                    [
                        self.install_button,
                        self.post_install_button,
                        self.update_button,
                        self.remove_button,
                        self.clear_output_button,
                    ]
                ),
                self.message_container,
                self.output_container,
            ],
        )

    def reconcile_state(self) -> None:
        """Refresh package compatibility and every widget derived from it."""
        self.installed_version, self.plugin_compatible, self.requirement_error = (
            get_plugin_version_info(
                self.data.package, self.data.pip or self.data.package
            )
        )
        self.app_compatible = is_version_compatible(self.data.requires_aiidalab_qe)
        if self.version_warning is None:
            return

        warning = ""
        if not self.app_compatible:
            warning += (
                '<div class="alert alert-danger" role="alert">'
                f"This plugin requires aiidalab_qe {self.data.requires_aiidalab_qe}, "
                f"but you have {INSTALLED_AIIDA_QE_VERSION}.</div>"
            )
        if self.requirement_error:
            warning += (
                '<div class="alert alert-danger" role="alert">'
                "Could not evaluate plugin requirement/version: "
                f"{self.requirement_error}</div>"
            )
        elif self.is_installed and not self.plugin_compatible:
            requirement = Requirement(self.data.pip or self.data.package)
            warning += (
                '<div class="alert alert-danger" role="alert">'
                f"Installed version {self.installed_version} does not satisfy "
                f"{requirement.name}{requirement.specifier}. Please update accordingly."
                "</div>"
            )
        self.version_warning.value = warning

        if self.install_button is not None:
            self.install_button.disabled = (
                self.is_installed
                or not self.app_compatible
                or bool(self.requirement_error)
            )
        if self.update_button is not None:
            self.update_button.disabled = (
                not self.is_installed
                or self.plugin_compatible
                or not self.app_compatible
            )
        if self.remove_button is not None:
            self.remove_button.disabled = not self.is_installed
        self._sync_title(warning)

    def _details_html(self) -> str:
        status_value = self.data.status.strip().lower()
        badge_color = COLOR_MAP.get(status_value, "#666666")
        display_text = status_value.capitalize() if status_value else "N/A"
        badge_html = (
            f'<span style="background-color: {badge_color}; color: #FFFFFF; '
            f'border-radius: 4px; padding: 2px 6px;">{display_text}</span>'
        )
        details = (
            f"<b>Package:</b> {self.data.package}<br>"
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

    def _sync_title(self, warning: str | None = None) -> None:
        if (
            self.accordion is None
            or self.index is None
            or self.index >= len(self.accordion.children)
        ):
            return
        if warning is None and self.version_warning is not None:
            warning = self.version_warning.value
        status = (
            "⚠️" if self.is_installed and warning else "✅" if self.is_installed else ""
        )
        self.accordion.set_title(self.index, f"{self.data.title} {status}".rstrip())

    def _on_install(self, _button: ipw.Button) -> None:
        message = self.message_container
        output = self.output_container
        output.layout.display = "block"
        message.layout.display = "block"
        message.value = ""
        output.value = ""
        self._append_message(f"Installing {self.data.package}...")
        install_source = self.data.pip or f"git+{self.data.github}"
        installed = self._execute_command(
            [sys.executable, "-m", "pip", "install", install_source, "--user"],
        )
        if not installed:
            self._append_message(
                "Installation did not complete. Review the command output for details.",
                color="#FF0000",
            )
            self.reconcile_state()
            return

        if self.data.post_install and not self._run_post_install(clear_output=False):
            self.reconcile_state()
            return

        self._append_message("Testing plugin loading...", color="#008000")
        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "aiidalab_qe",
                    "test-plugin",
                    self.data.package,
                ],
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError:
            LOGGER.exception(
                "Could not run plugin loading test for %s", self.data.package
            )
            self._append_message(
                "The plugin test could not be started. Details were logged for debugging.",
                color="#FF0000",
            )
            self.reconcile_state()
            return

        if result.stdout:
            self._append_output(result.stdout)
        if result.stderr:
            self._append_output(result.stderr)
        if result.returncode == 0:
            self._append_message("Plugin test passed.", color="#008000")
            self._append_message("Plugin installed successfully.", color="#008000")
            self._restart_daemon()
        else:
            LOGGER.error(
                "Plugin test failed for %s (exit code %s): %s%s",
                self.data.package,
                result.returncode,
                result.stdout,
                result.stderr,
            )
            self._append_message(
                f"The plugin test for {self.data.package} did not pass. The package will be "
                "removed to prevent use in an incomplete state. Details were logged.",
                color="#FF0000",
            )
            self._remove_package(clear_output=False)
        self.reconcile_state()

    def _on_update(self, _button: ipw.Button) -> None:
        output = self.output_container
        message = self.message_container
        output.layout.display = "block"
        message.layout.display = "block"
        message.value = ""
        self._append_message(f"Updating {self.data.package}...")
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
        if result:
            self._restart_daemon()
        self.reconcile_state()
        if result and self.is_installed and self.plugin_compatible:
            self._append_message(
                f"Updated {self.data.package} to {self.installed_version}.",
                color="#008000",
            )
        elif result:
            reason = self.requirement_error or (
                f"Installed version {self.installed_version or 'unknown'} does not meet "
                "the required version."
            )
            self._append_message(reason, color="#FF0000")

    def _on_remove(self, _button: ipw.Button) -> None:
        self.message_container.layout.display = "block"
        self.output_container.layout.display = "block"
        self._remove_package()

    def _remove_package(self, clear_output: bool = True) -> bool:
        self._append_message(f"Removing {self.data.package}...")
        requirement = Requirement(self.data.pip or self.data.package)
        result = self._execute_command(
            [sys.executable, "-m", "pip", "uninstall", "-y", requirement.name],
            clear_output=clear_output,
        )
        if result:
            self._append_message(
                f"{self.data.package} removed successfully.",
                color="#008000",
            )
            self._restart_daemon()
        self.reconcile_state()
        return result

    def _on_post_install(self, _button: ipw.Button) -> None:
        self._run_post_install()
        self.reconcile_state()

    def _run_post_install(self, clear_output: bool = True) -> bool:
        output = self.output_container
        message = self.message_container
        output.layout.display = "block"
        message.layout.display = "block"
        if clear_output:
            output.value = ""
            message.value = ""
        self._append_message(f"Running post-install for {self.data.package}...")
        command = [
            sys.executable,
            "-m",
            self.data.package.replace("-", "_"),
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
            LOGGER.exception("Could not start post-install for %s", self.data.package)
            self._append_message(
                "Post-install could not be started. Details were logged for debugging.",
                color="#FF0000",
            )
            return False
        if result.stdout:
            self._append_output(result.stdout)
        if result.returncode == 0:
            self._append_message(
                f"Post-install completed for {self.data.package}.",
                color="#008000",
            )
            return True
        LOGGER.error(
            "Post-install for %s failed with exit code %s: %s",
            self.data.package,
            result.returncode,
            result.stdout,
        )
        self._append_message(
            "Post-install did not complete. The package is present, but setup may be "
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
    def _restart_daemon() -> None:
        subprocess.run(["verdi", "daemon", "restart"], capture_output=True, check=False)


class PluginManager:
    """
    A manager class that reads a plugin configuration file (YAML),
    and creates an interactive Accordion UI for installing/uninstalling
    those plugins in a Jupyter environment.
    """

    def __init__(self, config_source: str = DEFAULT_PLUGIN_CONFIG_SOURCE):
        """
        Initialize the PluginManager with a path to a YAML config or a URL.

        :param config_source: Either a local YAML file path or a URL to a remote YAML file.
        """
        self.config_source = config_source
        self.config_error = None
        try:
            self.plugins = {
                name: QeAppPluginData.from_mapping(name, data)
                for name, data in self._load_config().items()
            }
        except (
            OSError,
            requests.RequestException,
            ValueError,
            yaml.YAMLError,
        ) as error:
            LOGGER.exception("Could not load plugin registry from %s", config_source)
            self.config_error = str(error)
            self.plugins = {}

        self.accordion = ipw.Accordion()
        self.rows: dict[str, QeAppPluginRow] = {}

    @property
    def data(self) -> dict:
        """Return normalized registry entries for existing consumers."""
        return {
            name: plugin_data.as_mapping() for name, plugin_data in self.plugins.items()
        }

    def display_ui(self) -> None:
        """
        Display the Accordion UI in a Jupyter notebook.
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
        display(self.accordion)

    def _load_config(self) -> dict:
        """Load YAML from the configured local path or HTTP(S) URL."""
        try:
            if urlparse(self.config_source).scheme in {"http", "https"}:
                response = requests.get(self.config_source, timeout=10)
                response.raise_for_status()
                content = response.text
            else:
                content = Path(self.config_source).read_text(encoding="utf-8")
            data = yaml.safe_load(content)
        except requests.RequestException as error:
            raise ValueError(f"Could not fetch plugin registry: {error}") from error
        except (OSError, yaml.YAMLError) as error:
            raise ValueError(f"Could not read plugin registry: {error}") from error

        if data is None:
            return {}
        if not isinstance(data, Mapping):
            raise TypeError("plugin registry must be a mapping of names to entries")
        return data

    def _build_ui(self) -> None:
        """Build the Accordion UI from the validated plugin records."""
        self.rows = {
            name: QeAppPluginRow(plugin_data)
            for name, plugin_data in self.plugins.items()
        }
        panels = [
            row.build_panel(self.accordion, index)
            for index, row in enumerate(self.rows.values())
        ]
        self.accordion.children = panels
        for row in self.rows.values():
            row._sync_title()
