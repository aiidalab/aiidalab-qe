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
from importlib import metadata
from threading import Thread

import ipywidgets as ipw
import requests
import yaml
from IPython.display import display
from packaging.requirements import InvalidRequirement, Requirement
from packaging.specifiers import SpecifierSet
from packaging.version import InvalidVersion, Version, parse

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
        installed_version = parse(INSTALLED_AIIDA_QE_VERSION)

        # If the installed version is a pre-release, allow it if it matches the specifier
        if installed_version.is_prerelease:
            return installed_version in specifier or specifier.contains(
                installed_version, prereleases=True
            )
        else:
            return installed_version in specifier

    except Exception as e:
        print(f"Warning: Could not parse version requirement '{required_version}': {e}")
        return False


def fetch_plugins_from_github(url: str) -> dict:
    """
    Fetch the plugins.yaml file from GitHub and parse it into a dictionary.

    :param url: URL of the plugins.yaml file.
    :return: Parsed dictionary of plugins data.
    """
    import requests

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise error if request fails
        return yaml.safe_load(response.text)
    except requests.RequestException as e:
        print(f"⚠️ Warning: Failed to fetch plugins.yaml from GitHub: {e}")
        return {}


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
    plugin_name: str, pip_requirement: str | None
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
        version, prereleases=True
    )
    return installed_version, compatible, None


def stream_output(process: subprocess.Popen, output_widget: ipw.HTML) -> None:
    """
    Reads output from the process and forwards it to the output widget.
    """
    while True:
        output = process.stdout.readline()
        if process.poll() is not None and output == "":
            break
        if output:
            _append_output(output_widget, output)


def _append_output(output_widget: ipw.HTML, text: str) -> None:
    """Append escaped text while preserving terminal line breaks."""
    output_widget.value += (
        '<div style="background-color: #3B3B3B; color: #FFFFFF; '
        'white-space: pre-wrap; margin: 0; padding: 0 4px;">'
        + html.escape(text)
        + "</div>"
    )


def _append_message(
    message_widget: ipw.HTML, text: str, color: str = "#000000"
) -> None:
    """Append an escaped high-level status message."""
    message_widget.value += (
        f'<div style="color: {color}; padding: 0 4px;">{html.escape(text)}</div>'
    )


def _clear_plugin_output(message_widget: ipw.HTML, output_widget: ipw.HTML) -> None:
    """Clear and hide the high-level and command output widgets."""
    message_widget.value = ""
    output_widget.value = ""
    message_widget.layout.display = "none"
    output_widget.layout.display = "none"


def _sync_clear_output_button(
    message_widget: ipw.HTML, output_widget: ipw.HTML, button: ipw.Button
) -> None:
    """Enable clearing only when either output widget contains text."""
    button.disabled = not (message_widget.value or output_widget.value)


def _set_accordion_status(accordion: ipw.Accordion, index: int, status: str) -> None:
    """Replace the status marker in an accordion title."""
    title = accordion.get_title(index)
    for marker in ("✅", "⚠️", "☐"):
        suffix = f" {marker}"
        if title.endswith(suffix):
            title = title[: -len(suffix)]
            break
    accordion.set_title(index, f"{title} {status}".rstrip())


def execute_command_with_output(
    command: list,
    output_widget: ipw.HTML,
    install_btn: ipw.Button,
    remove_btn: ipw.Button,
    action: str = "install",
    clear_output: bool = True,
) -> bool:
    """
    Execute a shell command and stream its output to the provided widget.
    """
    if clear_output:
        output_widget.value = ""
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
        _append_output(
            output_widget,
            "The command could not be started. Details were logged for debugging.\n",
        )
        return False

    thread = Thread(target=stream_output, args=(process, output_widget))
    thread.start()
    thread.join()

    if process.returncode == 0 and action == "install":
        install_btn.disabled = True
        remove_btn.disabled = False
        return True
    elif process.returncode == 0 and action == "update":
        install_btn.disabled = True
        remove_btn.disabled = False
        return True
    elif process.returncode == 0 and action == "remove":
        install_btn.disabled = False
        remove_btn.disabled = True
        return True
    else:
        LOGGER.error(
            "Plugin command %r failed with exit code %s", command, process.returncode
        )
        _append_output(
            output_widget,
            "\nThe command did not complete successfully. Review the output above; "
            "details were logged for debugging.\n",
        )
        return False


def remove_package(
    package_name: str,
    output_container: ipw.HTML,
    message_container: ipw.HTML,
    install_btn: ipw.Button,
    remove_btn: ipw.Button,
    accordion: ipw.Accordion,
    index: int,
    clear_output: bool = True,
) -> None:
    """
    Remove a plugin package via pip uninstall.
    """
    output_container.layout.display = "block"
    message_container.layout.display = "block"
    _append_message(message_container, f"Removing {package_name}...")
    normalized_name = package_name.replace("-", "_")
    command = ["pip", "uninstall", "-y", normalized_name]
    result = execute_command_with_output(
        command,
        output_container,
        install_btn,
        remove_btn,
        action="remove",
        clear_output=clear_output,
    )

    if result:
        _append_message(
            message_container,
            f"{package_name} removed successfully.",
            color="#008000",
        )
        _set_accordion_status(accordion, index, "")

        # Attempt to restart AiiDA daemon
        command = ["verdi", "daemon", "restart"]
        subprocess.run(command, capture_output=True, check=False)


def install_package(
    package_name: str,
    pip_install: str,
    github: str,
    post_install: str,
    output_container: ipw.HTML,
    message_container: ipw.HTML,
    install_btn: ipw.Button,
    remove_btn: ipw.Button,
    accordion: ipw.Accordion,
    index: int,
) -> None:
    """
    Install a plugin package from pip or GitHub, then optionally run a post-install command
    and test the plugin.
    """
    output_container.layout.display = "block"
    message_container.layout.display = "block"
    message_container.value = ""
    output_container.value = ""
    _append_message(message_container, f"Installing {package_name}...")

    if pip_install:
        command = ["pip", "install", pip_install, "--user"]
    else:
        command = ["pip", "install", "git+" + github, "--user"]

    install_result = execute_command_with_output(
        command, output_container, install_btn, remove_btn
    )
    if not install_result:
        _append_message(
            message_container,
            "Installation did not complete. Review the command output for details.",
            color="#FF0000",
        )
        return

    if post_install:
        _append_message(
            message_container, "Post-install step in progress...", color="#008000"
        )
        if not run_post_install(
            package_name,
            post_install,
            output_container,
            message_container,
            clear_output=False,
        ):
            _set_accordion_status(accordion, index, "⚠️")
            return

    _append_message(message_container, "Testing plugin loading...", color="#008000")
    cmd = [sys.executable, "-m", "aiidalab_qe", "test-plugin", package_name]
    try:
        test_result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        LOGGER.exception("Could not run plugin loading test for %s", package_name)
        _append_message(
            message_container,
            "The plugin test could not be started. Details were logged for debugging.",
            color="#FF0000",
        )
        _set_accordion_status(accordion, index, "⚠️")
        return

    if test_result.stdout:
        _append_output(output_container, test_result.stdout)
    if test_result.stderr:
        _append_output(output_container, test_result.stderr)

    if test_result.returncode == 0:
        _append_message(message_container, "Plugin test passed.", color="#008000")
        _append_message(
            message_container,
            "Plugin installed successfully.",
            color="#008000",
        )
        _set_accordion_status(accordion, index, "✅")

        # Restart daemon
        daemon_cmd = ["verdi", "daemon", "restart"]
        subprocess.run(daemon_cmd, capture_output=True, check=False)

    else:
        LOGGER.error(
            "Plugin test failed for %s (exit code %s): %s%s",
            package_name,
            test_result.returncode,
            test_result.stdout,
            test_result.stderr,
        )
        _append_message(
            message_container,
            f"The plugin test for {package_name} did not pass. The package will be "
            "removed to prevent use in an incomplete state. Details were logged.",
            color="#FF0000",
        )
        remove_package(
            package_name,
            output_container,
            message_container,
            install_btn,
            remove_btn,
            accordion,
            index,
            clear_output=False,
        )


def update_package(
    package_name: str,
    pip_install: str,
    github: str,
    output_container: ipw.HTML,
    message_container: ipw.HTML,
    install_btn: ipw.Button,
    update_btn: ipw.Button,
    remove_btn: ipw.Button,
    version_warning: ipw.HTML,
    accordion: ipw.Accordion | None = None,
    index: int | None = None,
) -> None:
    """Upgrade a plugin and clear its warning if it meets the registry requirement."""
    output_container.layout.display = "block"
    message_container.layout.display = "block"
    message_container.value = ""
    _append_message(message_container, f"Updating {package_name}...")
    requirement = pip_install or "git+" + github
    result = execute_command_with_output(
        ["pip", "install", "--upgrade", requirement, "--user"],
        output_container,
        install_btn,
        remove_btn,
        action="update",
    )

    if result:
        installed_version, compatible, error = get_plugin_version_info(
            package_name, pip_install
        )
        if compatible and installed_version is not None:
            version_warning.value = ""
            update_btn.disabled = True
            if accordion is not None and index is not None:
                _set_accordion_status(accordion, index, "✅")
            _append_message(
                message_container,
                f"Updated {package_name} to {installed_version}.",
                color="#008000",
            )
        elif error:
            _append_message(
                message_container,
                f"Could not verify the installed version: {error}",
                color="#FF0000",
            )
        else:
            _append_message(
                message_container,
                f"Installed version {installed_version or 'unknown'} does not meet "
                "the required version.",
                color="#FF0000",
            )


def run_post_install(
    package_name: str,
    post_install: str,
    output_container: ipw.HTML,
    message_container: ipw.HTML,
    clear_output: bool = True,
) -> bool:
    """Run only the configured post-install command for an installed plugin."""
    output_container.layout.display = "block"
    message_container.layout.display = "block"
    if clear_output:
        output_container.value = ""
        message_container.value = ""
    _append_message(message_container, f"Running post-install for {package_name}...")
    command = [sys.executable, "-m", package_name.replace("-", "_"), post_install]
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    except OSError:
        LOGGER.exception("Could not start post-install for %s", package_name)
        _append_message(
            message_container,
            "Post-install could not be started. Details were logged for debugging.",
            color="#FF0000",
        )
        return False
    if result.stdout:
        _append_output(output_container, result.stdout)
    if result.returncode == 0:
        _append_message(
            message_container,
            f"Post-install completed for {package_name}.",
            color="#008000",
        )
        return True
    else:
        LOGGER.error(
            "Post-install for %s failed with exit code %s: %s",
            package_name,
            result.returncode,
            result.stdout,
        )
        _append_message(
            message_container,
            "Post-install did not complete. The package is present, but setup may be "
            "incomplete. Review the command output; details were logged for debugging.",
            color="#FF0000",
        )
        return False


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
        self.config_source = "/home/jovyan/apps/quantum-espresso/plugins.yaml"
        self.data = self._load_config()
        self.accordion = ipw.Accordion()

    def _load_config(self) -> dict:
        """
        Load and parse the YAML plugin configuration from a local file or a URL.
        """
        if self.config_source.startswith("http"):  # Detect if it's a URL
            return self._fetch_remote_config(self.config_source)
        else:  # Assume it's a local file
            return self._fetch_local_config(self.config_source)

    def _fetch_remote_config(self, url: str) -> dict:
        """
        Fetch the plugins.yaml file from GitHub and parse it into a dictionary.

        :param url: URL of the plugins.yaml file.
        :return: Parsed dictionary of plugins data.
        """
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()  # Raise an error if request fails
            return yaml.safe_load(response.text)
        except requests.RequestException as e:
            print(f"⚠️ Warning: Failed to fetch plugins.yaml from GitHub: {e}")
            return {}

    def _fetch_local_config(self, file_path: str) -> dict:
        """
        Load plugins.yaml from a local file.

        :param file_path: Path to the YAML file.
        :return: Parsed dictionary of plugins data.
        """
        try:
            with open(file_path) as file:
                return yaml.safe_load(file)
        except (FileNotFoundError, yaml.YAMLError) as e:
            print(f"⚠️ Warning: Failed to load local plugins.yaml: {e}")
            return {}

    def _build_ui(self) -> None:
        """
        Build the Accordion UI based on the loaded plugin data.
        """
        for i, (plugin_name, plugin_data) in enumerate(self.data.items()):
            package_name = plugin_data.get("package", plugin_name)
            pip_data = plugin_data.get("pip")
            installed_version, plugin_compatible, requirement_error = (
                get_plugin_version_info(package_name, pip_data or package_name)
            )
            installed = installed_version is not None
            required_version = plugin_data.get("requires_aiidalab_qe", None)
            app_compatible = is_version_compatible(required_version)

            # Create message and output containers
            message_container = ipw.HTML(
                value="",
                layout=ipw.Layout(
                    max_height="250px",
                    overflow="auto",
                    border="1px solid #9e9e9e",
                    display="none",
                ),
            )
            output_container = ipw.HTML(
                value="",
                layout=ipw.Layout(
                    max_height="250px",
                    overflow="auto",
                    border="1px solid #9e9e9e",
                    display="none",
                ),
            )

            status_value = plugin_data.get("status", "").strip().lower()
            badge_color = COLOR_MAP.get(status_value, "#666666")
            display_text = status_value.capitalize() if status_value else "N/A"

            badge_html = f"""
                <span style="background-color: {badge_color}; color: #FFFFFF;
                            border-radius: 4px; padding: 2px 6px;">
                    {display_text}
                </span>
            """

            # Display host app and plugin package version requirements.
            version_message = ""
            if not app_compatible:
                version_message = f"""
                <div class="alert alert-danger" role="alert">
                    ⚠️ This plugin requires aiidalab_qe {required_version}, but you have {INSTALLED_AIIDA_QE_VERSION}.
                </div>
                """
            if requirement_error:
                version_message += f"""
                <div class="alert alert-danger" role="alert">
                    ⚠️ Could not evaluate plugin requirement/version: {requirement_error}
                </div>
                """
            elif installed and not plugin_compatible:
                requirement = Requirement(pip_data or package_name)
                version_message += f"""
                <div class="alert alert-danger" role="alert">
                    ⚠️ Installed version {installed_version} does not satisfy
                    {requirement.name}{requirement.specifier}. Please update accordingly.
                </div>
                """

            # Build plugin description details
            details = f"""
                <b>Package:</b> {package_name}<br>
                <b>Author:</b> {plugin_data.get("author", "N/A")}<br>
                <b>Description:</b> {plugin_data.get("description", "No description available")}<br>
                <b>Status:</b> {badge_html}<br>
            """

            if "documentation" in plugin_data:
                details += f"<b>Documentation:</b> <a href='{plugin_data['documentation']}' target='_blank'>Visit</a><br>"
            if "github" in plugin_data:
                details += f"<b>Github:</b> <a href='{plugin_data.get('github')}' target='_blank'>Visit</a>"

            # Create install/remove buttons
            install_btn = ipw.Button(
                description="Install",
                button_style="success",
                disabled=installed or not app_compatible or bool(requirement_error),
            )
            update_btn = ipw.Button(
                description="Update",
                button_style="warning",
                disabled=not installed or plugin_compatible or not app_compatible,
            )
            remove_btn = ipw.Button(
                description="Remove",
                button_style="danger",
                disabled=not installed,
            )
            clear_output_btn = ipw.Button(
                description="Clear output",
                icon="trash-o",
                tooltip="Clear plugin messages and command output",
                disabled=True,
            )
            post_install_btn = ipw.Button(
                description="Run post-install only",
                button_style="",
            )
            post_install_btn.layout.display = (
                "" if plugin_data.get("post_install") else "none"
            )
            version_warning = ipw.HTML(value=version_message)

            # Attach callbacks
            github_data = plugin_data.get("github", "")
            post_install_data = plugin_data.get("post_install", None)
            install_btn.on_click(
                lambda _btn, pn=package_name, pip=pip_data, gh=github_data, post=post_install_data, oc=output_container, mc=message_container, ib=install_btn, rb=remove_btn, ac=self.accordion, idx=i: (
                    install_package(pn, pip, gh, post, oc, mc, ib, rb, ac, idx)
                )
            )
            remove_btn.on_click(
                lambda _btn, pn=package_name, oc=output_container, mc=message_container, ib=install_btn, rb=remove_btn, ac=self.accordion, idx=i: (
                    remove_package(pn, oc, mc, ib, rb, ac, idx)
                )
            )
            update_btn.on_click(
                lambda _btn, pn=package_name, pip=pip_data, gh=github_data, oc=output_container, mc=message_container, ib=install_btn, ub=update_btn, rb=remove_btn, warning=version_warning, ac=self.accordion, idx=i: (
                    update_package(pn, pip, gh, oc, mc, ib, ub, rb, warning, ac, idx)
                )
            )
            post_install_btn.on_click(
                lambda _btn, pn=package_name, post=post_install_data, oc=output_container, mc=message_container: (
                    run_post_install(pn, post, oc, mc)
                )
            )
            clear_output_btn.on_click(
                lambda _btn, oc=output_container, mc=message_container: (
                    _clear_plugin_output(mc, oc)
                )
            )
            for output_widget in (message_container, output_container):
                output_widget.observe(
                    lambda _change, mc=message_container, oc=output_container, btn=clear_output_btn: (
                        _sync_clear_output_button(mc, oc, btn)
                    ),
                    names="value",
                )

            # Create layout for each plugin
            box = ipw.VBox(
                [
                    ipw.HTML(details),
                    version_warning,
                    ipw.HBox(
                        [
                            install_btn,
                            update_btn,
                            post_install_btn,
                            remove_btn,
                            clear_output_btn,
                        ]
                    ),
                    message_container,
                    output_container,
                ]
            )

            # Keep the title marker consistent with the compatibility warning.
            status_icon = (
                "⚠️" if installed and version_message else "✅" if installed else ""
            )
            title_with_icon = f"{plugin_data.get('title')} {status_icon}".rstrip()
            self.accordion.children = [*self.accordion.children, box]
            self.accordion.set_title(i, title_with_icon)

    def display_ui(self) -> None:
        """
        Display the Accordion UI in a Jupyter notebook.
        """
        self.accordion.children = ()
        self._build_ui()
        display(self.accordion)
