"""
plugin_manager.py

Module that provides a plugin management interface, allowing the user to
install and remove Python-based AiiDAlab plugins, with real-time streaming
of command output in a Jupyter environment.
"""

import importlib
import subprocess
import sys
from pathlib import Path
from threading import Thread

import ipywidgets as ipw
import yaml
from IPython.display import display

# Define badge colors based on status
COLOR_MAP = {
    "experimental": "#FF8C00",  # 🟠 Orange - Early development
    "beta": "#FFA500",  # 🟡 Darker Orange - Testing phase
    "stable": "#1E90FF",  # 🔵 Dodger Blue - Well-tested, feature-complete
    "production": "#008000",  # 🟢 Green - Fully stable and recommended
    "deprecated": "#FF0000",  # 🔴 Red - No longer maintained
    "archived": "#808080",  # ⚪ Grey - Retained for reference, no updates
}


def is_package_installed(package_name: str) -> bool:
    """
    Check if a given Python package is already installed.

    :param package_name: Name of the package to check.
    :return: True if installed, False otherwise.
    """
    package_name = package_name.replace("-", "_")
    try:
        importlib.import_module(package_name)
    except ImportError:
        return False
    else:
        return True


def stream_output(process: subprocess.Popen, output_widget: ipw.HTML) -> None:
    """
    Reads output from the process and forwards it to the output widget.
    """
    while True:
        output = process.stdout.readline()
        if process.poll() is not None and output == "":
            break
        if output:
            output_widget.value += f"""<div style="background-color: #3B3B3B; color: #FFFFFF;">{output}</div>"""


def execute_command_with_output(
    command: list,
    output_widget: ipw.HTML,
    install_btn: ipw.Button,
    remove_btn: ipw.Button,
    action: str = "install",
) -> bool:
    """
    Execute a shell command and stream its output to the provided widget.
    """
    output_widget.value = ""  # Clear the widget before running the command
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )

    thread = Thread(target=stream_output, args=(process, output_widget))
    thread.start()
    thread.join()

    if process.returncode == 0 and action == "install":
        install_btn.disabled = True
        remove_btn.disabled = False
        return True
    elif process.returncode == 0 and action == "remove":
        install_btn.disabled = False
        remove_btn.disabled = True
        return True
    else:
        output_widget.value += """<div style="background-color: #3B3B3B; color: #FF0000;">Command failed.</div>"""
        return False


def remove_package(
    package_name: str,
    output_container: ipw.HTML,
    message_container: ipw.HTML,
    install_btn: ipw.Button,
    remove_btn: ipw.Button,
    accordion: ipw.Accordion,
    index: int,
) -> None:
    """
    Remove a plugin package via pip uninstall.
    """
    # Show the containers when user starts the remove action.
    message_container.layout.display = "block"
    output_container.layout.display = "block"

    message_container.value += (
        f"""<div style="color: #FF0000;">Removing {package_name}...</div>"""
    )
    normalized_name = package_name.replace("-", "_")
    command = ["pip", "uninstall", "-y", normalized_name]
    result = execute_command_with_output(
        command, output_container, install_btn, remove_btn, action="remove"
    )

    if result:
        message_container.value += f"""<div style="color: #008000;">{package_name} removed successfully.</div>"""
        # Remove the checkmark from the Accordion title
        title = accordion.get_title(index)
        if title.endswith("✅"):
            new_title = title[:-2] + "☐"
            accordion.set_title(index, new_title)

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
    # Show the containers when user starts the install action.
    message_container.layout.display = "block"
    output_container.layout.display = "block"

    if pip_install:
        command = ["pip", "install", pip_install, "--user"]
    else:
        command = ["pip", "install", "git+" + github, "--user"]

    message_container.value = (
        f"""<div style="color: #000000;">Installing {package_name}...</div>"""
    )
    install_result = execute_command_with_output(
        command, output_container, install_btn, remove_btn
    )

    # If post_install is defined, attempt to run that command
    if install_result and post_install:
        message_container.value += (
            """<div style="color: #008000;">Post-install step in progress...</div>"""
        )
        cmd = [sys.executable, "-m", package_name.replace("-", "_"), post_install]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            message_container.value += f"""
                <div style="color: #FF0000;">
                    Post-install command failed: {result.stderr}
                </div>
            """

    # Test the plugin functionality if install was successful
    if install_result:
        message_container.value += (
            """<div style="color: #008000;">Testing plugin loading...</div>"""
        )
        cmd = [sys.executable, "-m", "aiidalab_qe", "test-plugin", package_name]
        test_result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if test_result.returncode == 0:
            message_container.value += (
                """<div style="color: #008000;">Plugin test passed.</div>"""
            )
            message_container.value += (
                """<div style="color: #008000;">Plugin installed successfully.</div>"""
            )

            # Update Accordion title with a checkmark
            title = accordion.get_title(index)
            if title.endswith("☐"):
                new_title = title[:-1] + "✅"
                accordion.set_title(index, new_title)

            # Restart daemon
            daemon_cmd = ["verdi", "daemon", "restart"]
            subprocess.run(daemon_cmd, capture_output=True, check=False)

        else:
            message_container.value += f"""
                <div style="color: #FF0000;">
                    The plugin '{package_name}' was installed but failed functionality test:
                    {test_result.stderr}.
                </div>
            """
            message_container.value += """
                <div style="color: #FF0000;">
                    This may be due to compatibility issues with the current AiiDAlab QEApp.
                    Please contact the plugin author for assistance.
                </div>
            """
            message_container.value += """
                <div style="color: #FF0000;">
                    The plugin will now be uninstalled to prevent issues.
                </div>
            """
            remove_package(
                package_name,
                output_container,
                message_container,
                install_btn,
                remove_btn,
                accordion,
                index,
            )


class PluginManager:
    """
    A manager class that reads a plugin configuration file (YAML),
    and creates an interactive Accordion UI for installing/uninstalling
    those plugins in a Jupyter environment.
    """

    def __init__(self, config_path: str = "plugins.yaml"):
        """
        Initialize the PluginManager with a path to a YAML config.
        """
        self.config_path = Path(config_path).resolve()
        self.data = self._load_config()
        self.accordion = ipw.Accordion()
        self._build_ui()

    def _load_config(self) -> dict:
        """
        Load and parse the YAML plugin configuration.
        """
        with self.config_path.open("r") as file:
            return yaml.safe_load(file)

    def _build_ui(self) -> None:
        """
        Build the Accordion UI based on the loaded plugin data.
        """
        for i, (plugin_name, plugin_data) in enumerate(self.data.items()):
            installed = is_package_installed(plugin_name)

            # Create message and output containers
            output_container = ipw.HTML(
                value="""
                <div style="background-color: #3B3B3B; color: #FFFFFF; height: 100%; overflow: auto;">
                </div>
                """,
                layout=ipw.Layout(
                    max_height="250px",
                    overflow="auto",
                    border="2px solid #CCCCCC",
                    display="none",  # hide by default
                ),
            )
            message_container = ipw.HTML(
                value="""
                <div style="color: #000000; height: 100%; overflow: auto;">
                </div>
                """,
                layout=ipw.Layout(
                    max_height="250px",
                    overflow="auto",
                    border="2px solid #CCCCCC",
                    display="none",  # hide by default
                ),
            )

            status_value = plugin_data.get("status", "").strip().lower()
            # Default color if no status or unknown status
            badge_color = COLOR_MAP.get(status_value, "#666666")
            if status_value:
                # Capitalize first letter for display
                display_text = status_value.capitalize()
            else:
                display_text = "N/A"

            badge_html = f"""
                <span style="background-color: {badge_color}; color: #FFFFFF;
                             border-radius: 4px; padding: 2px 6px;">
                    {display_text}
                </span>
            """

            # Build plugin description details
            details = f"Author: {plugin_data.get('author', 'N/A')}<br>"
            details += f"Description: {plugin_data.get('description', 'No description available')}<br>"
            details += f"Status: {badge_html}<br>"  # <<--- Add badge here

            if "documentation" in plugin_data:
                details += f"Documentation: <a href='{plugin_data['documentation']}' target='_blank'>Visit</a><br>"
            if "github" in plugin_data:
                details += f"Github: <a href='{plugin_data.get('github')}' target='_blank'>Visit</a>"

            # Create install/remove buttons
            install_btn = ipw.Button(
                description="Install", button_style="success", disabled=installed
            )
            remove_btn = ipw.Button(
                description="Remove", button_style="danger", disabled=not installed
            )

            # Attach callbacks
            pip_data = plugin_data.get("pip", None)
            github_data = plugin_data.get("github", "")
            post_install_data = plugin_data.get("post_install", None)
            install_btn.on_click(
                lambda _btn,
                pn=plugin_name,
                pip=pip_data,
                gh=github_data,
                post=post_install_data,
                oc=output_container,
                mc=message_container,
                ib=install_btn,
                rb=remove_btn,
                ac=self.accordion,
                idx=i: install_package(pn, pip, gh, post, oc, mc, ib, rb, ac, idx)
            )
            remove_btn.on_click(
                lambda _btn,
                pn=plugin_name,
                oc=output_container,
                mc=message_container,
                ib=install_btn,
                rb=remove_btn,
                ac=self.accordion,
                idx=i: remove_package(pn, oc, mc, ib, rb, ac, idx)
            )

            # Create layout for each plugin
            box = ipw.VBox(
                [
                    ipw.HTML(details),
                    ipw.HBox([install_btn, remove_btn]),
                    message_container,
                    output_container,
                ]
            )

            # Define the title with checkmark if installed
            title_with_icon = f"{plugin_data.get('title')} {'✅' if installed else '☐'}"
            self.accordion.set_title(i, title_with_icon)
            self.accordion.children = [*self.accordion.children, box]

    def display_ui(self) -> None:
        """
        Display the Accordion UI in a Jupyter notebook.
        """
        display(self.accordion)
