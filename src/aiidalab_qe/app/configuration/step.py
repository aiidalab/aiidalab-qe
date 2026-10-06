"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

from importlib.metadata import distributions

import ipywidgets as ipw
from packaging.requirements import InvalidRequirement, Requirement
from packaging.utils import canonicalize_name

from aiidalab_qe.app.utils.plugin_manager import (
    DEFAULT_PLUGIN_CONFIG_SOURCE,
    PluginManager,
    get_plugin_version_info,
    is_package_installed,
    is_version_compatible,
)
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel, PanelModel
from aiidalab_qe.common.widgets import LinkButton
from aiidalab_qe.common.wizard import ConfirmableDependentWizardStep
from aiidalab_qe.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.plugins.state import get_activation_failure
from aiidalab_qe.plugins.utils import get_entry_items

from .advanced import (
    AdvancedConfigurationSettingsModel,
    AdvancedConfigurationSettingsPanel,
)
from .basic import BasicConfigurationSettingsModel, BasicConfigurationSettingsPanel
from .model import ConfigurationStepModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class ConfigurationStep(ConfirmableDependentWizardStep[ConfigurationStepModel]):
    _missing_message = "Missing input structure"

    def __init__(self, model: ConfigurationStepModel, **kwargs):
        super().__init__(
            model=model,
            confirm_kwargs={
                "tooltip": "Confirm the currently selected settings and go to the next step",
            },
            **kwargs,
        )

        workchain_model = BasicConfigurationSettingsModel()
        self.workchain_settings = BasicConfigurationSettingsPanel(model=workchain_model)
        self._model.add_model("workchain", workchain_model)

        advanced_model = AdvancedConfigurationSettingsModel()
        self.advanced_settings = AdvancedConfigurationSettingsPanel(
            model=advanced_model
        )
        self._model.add_model("advanced", advanced_model)

        # HACK due to spin orbit moving to basic settings (#984), we need to
        # sync the basic model's spin orbit when the advanced model's spin
        # orbit is preloaded
        ipw.dlink(
            (advanced_model, "spin_orbit"),
            (workchain_model, "spin_orbit"),
        )

        self._model.observe(
            self._on_input_structure_change,
            "structure_uuid",
        )
        self._model.observe(
            self._on_installed_properties_fetched,
            "installed_properties_fetched",
        )
        self._model.observe(
            self._on_available_properties_fetched,
            "available_properties_fetched",
        )

        self.settings = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        self.installed_properties_list = []
        self.incompatible_properties_list = []
        self.available_properties_list = []
        self.incompatible_plugin_data = {}
        self.properties = {}
        self.observed_property_models = set()

        self.entry_point_distributions = {
            entry_point.name: distribution.metadata["Name"]
            for distribution in distributions()
            for entry_point in distribution.entry_points
            if entry_point.group == "aiidalab_qe.properties"
            and distribution.metadata.get("Name")
        }

        self._fetch_available_properties()
        self._fetch_plugin_calculation_settings()

    def reset(self):
        self._model.reset()
        if self.rendered:
            self.sub_steps.selected_index = None
            self.tabs.selected_index = 0

    def _render(self):
        super()._render()

        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type_help = ipw.HTML()
        ipw.dlink(
            (self._model, "relax_type_help"),
            (self.relax_type_help, "value"),
        )
        self.relax_type = ipw.ToggleButtons()
        ipw.dlink(
            (self._model, "relax_type_options"),
            (self.relax_type, "options"),
        )
        ipw.link(
            (self._model, "relax_type"),
            (self.relax_type, "value"),
        )

        self.tabs = ipw.Tab(
            layout=ipw.Layout(min_height="250px"),
            selected_index=None,
        )
        self.tabs.observe(
            self._on_tab_change,
            "selected_index",
        )

        self.install_new_plugin_button = LinkButton(
            description="Plugin store",
            link="plugin_manager.ipynb",
            icon="puzzle-piece",  # More intuitive icon
            tooltip="Browse and install additional plugins from the Plugin Store",
            layout=ipw.Layout(margin="4px 0"),
        )

        self.installed_properties = ipw.VBox(children=self.installed_properties_list)
        self.incompatible_properties = ipw.HTML()
        self.incompatible_properties_section = ipw.VBox(
            children=[
                ipw.HTML("<hr>"),
                ipw.HTML("<h4>Incompatible</h4>"),
                ipw.HTML(
                    "<em>The following plugins require attention in the Plugin store:</em>"
                ),
                self.incompatible_properties,
            ],
            layout=ipw.Layout(
                display="" if self.incompatible_properties_list else "none"
            ),
        )

        self.available_properties = ipw.HTML("""
            <div class="loading" style="display: flex; align-items: center; font-size: unset;">
                Loading available properties
                <i class="fa fa-spinner fa-spin fa-2x fa-fw"></i>
            </div>
        """)

        self.sub_steps = ipw.Accordion(
            children=[
                ipw.VBox(
                    children=[
                        InAppGuide(identifier="properties-selection"),
                        ipw.HTML("<h4>Ready for use</h4>"),
                        ipw.HTML(
                            "<em>Select a property to add its settings panel in "
                            "step 2.2:</em>"
                        ),
                        self.installed_properties,
                        self.incompatible_properties_section,
                        ipw.HTML("<hr>"),
                        ipw.HTML("<h4>Available in store</h4>"),
                        ipw.HTML(
                            "<em>The following properties are available in the "
                            "plugin store:</em>"
                        ),
                        self.available_properties,
                        ipw.HTML("<hr>"),
                        ipw.HTML(
                            "Visit the plugin store to browse and install additional plugins, "
                            "or to resolve incompatible plugins.<br/>"
                            "<b>Note:</b> The app <b>must be reloaded</b> for changes to take "
                            "effect.",
                        ),
                        self.install_new_plugin_button,
                    ]
                ),
                ipw.VBox(
                    children=[
                        InAppGuide(identifier="calculation-settings"),
                        self.tabs,
                    ],
                ),
            ],
            layout=ipw.Layout(margin="10px 2px"),
            selected_index=None,
        )
        self.sub_steps.set_title(0, "Step 2.1: Select which properties to calculate")
        self.sub_steps.set_title(1, "Step 2.2: Customize calculation parameters")

        self.content.children = [
            InAppGuide(identifier="configuration-step"),
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Structure relaxation</h4>
                </div>
            """),
            self.relax_type_help,
            self.relax_type,
            self.sub_steps,
            self.confirm_box,
        ]

    def _post_render(self):
        super()._post_render()
        self._set_incompatible_properties()
        self._set_available_properties()
        self._update_tabs()

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: ConfigurationSettingsPanel = self.tabs.children[tab_index]  # type: ignore
        tab.render()

    def _on_input_structure_change(self, _):
        self._model.update()

    def _on_installed_properties_fetched(self, _):
        if not self.rendered:
            return
        self.installed_properties.children = self.installed_properties_list
        self._set_incompatible_properties()
        self._toggle_incompatible_properties()

    def _on_available_properties_fetched(self, _):
        if not self.rendered:
            return
        self._set_incompatible_properties()
        self._toggle_incompatible_properties()
        self._set_available_properties()

    def _toggle_incompatible_properties(self):
        self.incompatible_properties_section.layout.display = (
            "" if self.incompatible_properties_list else "none"
        )

    def _set_incompatible_properties(self):
        items = "".join(
            f"<li>{title}</li>" for title in self.incompatible_properties_list
        )
        self.incompatible_properties.value = (
            f"<ul style='margin: 0;'>{items}</ul>" if items else ""
        )

    def _set_available_properties(self):
        self.available_properties.value = f"""
            <ul style='margin: 0;'>
                {"".join(f"<li>{title}</li>" for title in self.available_properties_list)}
            </ul>
        """

    def _update_tabs(self):
        children = []
        titles = []
        for identifier, model in self._model.get_models():
            if model.include:
                settings = self.settings[identifier]
                titles.append(model.title)
                children.append(settings)
        if self.rendered:
            self.tabs.selected_index = None
            self.tabs.children = children
            for i, title in enumerate(titles):
                self.tabs.set_title(i, title)
            self.tabs.selected_index = 0

    def _fetch_plugin_calculation_settings(self):
        include_activation_failures = self._model.loaded_from_process
        outlines = get_entry_items(
            "aiidalab_qe.properties",
            "outline",
            include_activation_failures=include_activation_failures,
        )
        entries = get_entry_items(
            "aiidalab_qe.properties",
            "configuration",
            include_activation_failures=include_activation_failures,
        )

        self.incompatible_properties_list = [
            plugin_data["title"]
            for plugin_data in self.incompatible_plugin_data.values()
        ]

        for identifier, configuration in entries.items():
            for key in ("panel", "model"):
                if key not in configuration:
                    raise ValueError(f"Entry {identifier} is missing the '{key}' key")

            model: PanelModel = configuration["model"]()
            self._model.add_model(identifier, model)

            owner = self.entry_point_distributions.get(identifier)
            plugin_data = self.incompatible_plugin_data.get(
                canonicalize_name(owner) if owner else None
            )
            if plugin_data is not None:
                model.include = False
                self.settings[identifier] = configuration["panel"](model=model)
                continue

            outline = outlines[identifier]()
            info = ipw.HTML()
            ipw.link(
                (model, "include"),
                (outline.include, "value"),
            )

            if identifier == "bands":
                ipw.dlink(
                    (self._model, "structure_uuid"),
                    (outline.include, "disabled"),
                    lambda _: not self._model.has_pbc,
                )

            panel: ConfigurationSettingsPanel = configuration["panel"](model=model)

            def toggle_plugin(_, panel=panel):
                panel.refresh()
                self._update_tabs()

            model.observe(
                toggle_plugin,
                "include",
            )

            self.installed_properties_list.append(
                ipw.HBox(
                    children=[
                        outline,
                        info,
                    ]
                )
            )

            self.settings[identifier] = panel

        self._model.installed_properties_fetched = True

    def _fetch_available_properties(self, plugin_config_source=None):
        plugin_config_source = plugin_config_source or DEFAULT_PLUGIN_CONFIG_SOURCE
        plugin_manager = PluginManager(plugin_config_source)

        for plugin_name, plugin_data in plugin_manager.data.items():
            if (
                plugin_data.get("category", "calculation").lower() != "calculation"
            ):  # Ignore non-property plugins
                continue

            package_name = plugin_data.get("package") or plugin_name
            pip_requirement = plugin_data.get("pip") or package_name
            activation_error = (
                None
                if self._model.loaded_from_process
                else get_activation_failure(package_name)
            )

            if self._model.loaded_from_process:
                installed_version = None
                plugin_compatible = True
                requirement_error = None
                app_compatible = True
                is_installed = is_package_installed(package_name)
            else:
                (
                    installed_version,
                    plugin_compatible,
                    requirement_error,
                ) = get_plugin_version_info(package_name, pip_requirement)
                is_installed = installed_version is not None
                compatible_app_version = plugin_data.get("requires_aiidalab_qe")
                app_compatible = is_version_compatible(compatible_app_version)

            try:
                requirement = Requirement(pip_requirement)
                distribution_name = canonicalize_name(requirement.name)
                requirement_text = str(requirement)
            except InvalidRequirement:
                distribution_name = canonicalize_name(package_name)
                requirement_text = pip_requirement
                is_installed = is_package_installed(package_name)

            if is_installed and (
                not plugin_compatible
                or not app_compatible
                or requirement_error
                or activation_error
            ):
                self.incompatible_plugin_data[distribution_name] = {
                    **plugin_data,
                    "installed_version": installed_version,
                    "plugin_compatible": plugin_compatible,
                    "app_compatible": app_compatible,
                    "requirement_error": requirement_error,
                    "requirement": requirement_text,
                    "activation_error": activation_error,
                }

            if not is_installed:
                self.available_properties_list.append(plugin_data["title"])

        self._model.incompatible_properties = {
            identifier
            for identifier, owner in self.entry_point_distributions.items()
            if owner and canonicalize_name(owner) in self.incompatible_plugin_data
        }

        self._model.available_properties_fetched = True
