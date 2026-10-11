"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import asyncio
from importlib.metadata import distributions

import ipywidgets as ipw
from packaging.utils import canonicalize_name

from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel, PanelModel
from aiidalab_qe.common.widgets import LinkButton
from aiidalab_qe.common.wizard import ConfirmableDependentWizardStep
from aiidalab_qe.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.plugins.registry import (
    get_plugin_registry,
    get_plugin_status,
    get_plugin_statuses,
    is_plugin_installed,
)
from aiidalab_qe.plugins.utils import get_entries

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
            self._on_installed_plugins_fetched,
            "installed_plugins_fetched",
        )
        self._model.observe(
            self._on_available_plugins_fetched,
            "available_plugins_fetched",
        )

        self.settings = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        self.installed_plugins_list = []
        self.installed_plugin_rows = {}
        self.incompatible_plugins_list = []
        self.available_plugins_list = []
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

        self._fetch_available_plugins()
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

        self.installed_plugins = ipw.VBox(children=self.installed_plugins_list)
        self.incompatible_plugins = ipw.HTML()
        self.incompatible_plugins_section = ipw.VBox(
            children=[
                ipw.HTML("<hr>"),
                ipw.HTML("<h4>Incompatible</h4>"),
                ipw.HTML(
                    "<em>The following plugins require attention in the Plugin store:</em>"
                ),
                self.incompatible_plugins,
            ],
            layout=ipw.Layout(display="" if self.incompatible_plugins_list else "none"),
        )

        self.available_plugins = ipw.HTML("""
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
                        self.installed_plugins,
                        self.incompatible_plugins_section,
                        ipw.HTML("<hr>"),
                        ipw.HTML("<h4>Available in store</h4>"),
                        ipw.HTML(
                            "<em>The following plugins are available in the "
                            "plugin store:</em>"
                        ),
                        self.available_plugins,
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

        self._schedule_runtime_plugin_checks()

    def _post_render(self):
        super()._post_render()
        self._set_incompatible_plugins()
        self._set_available_plugins()
        self._update_tabs()

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: ConfigurationSettingsPanel = self.tabs.children[tab_index]  # type: ignore
        tab.render()

    def _on_input_structure_change(self, _):
        self._model.update()

    def _on_installed_plugins_fetched(self, _):
        if not self.rendered:
            return
        self.installed_plugins.children = self.installed_plugins_list
        self._set_incompatible_plugins()
        self._toggle_incompatible_plugins()

    def _on_available_plugins_fetched(self, _):
        if not self.rendered:
            return
        self._set_incompatible_plugins()
        self._toggle_incompatible_plugins()
        self._set_available_plugins()

    def _toggle_incompatible_plugins(self):
        self.incompatible_plugins_section.layout.display = (
            "" if self.incompatible_plugins_list else "none"
        )

    def _set_incompatible_plugins(self):
        items = "".join(f"<li>{title}</li>" for title in self.incompatible_plugins_list)
        self.incompatible_plugins.value = (
            f"<ul style='margin: 0;'>{items}</ul>" if items else ""
        )

    def _set_available_plugins(self):
        self.available_plugins.value = f"""
            <ul style='margin: 0;'>
                {"".join(f"<li>{title}</li>" for title in self.available_plugins_list)}
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
        loaded_entries = get_entries("aiidalab_qe.properties")
        outlines = {
            identifier: entry["outline"]
            for identifier, entry in loaded_entries.items()
            if entry.get("outline", False)
        }
        entries = {
            identifier: entry["configuration"]
            for identifier, entry in loaded_entries.items()
            if entry.get("configuration", False)
        }

        self.incompatible_plugins_list = [
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

            self.installed_plugin_rows[identifier] = outline
            self.installed_plugins_list.append(outline)

            self.settings[identifier] = panel

        self._model.installed_plugins_fetched = True

    def _fetch_available_plugins(
        self,
        *,
        statuses=None,
    ):
        registry = get_plugin_registry()
        self._plugin_registry = registry

        self.available_plugins_list = []
        self.incompatible_plugin_data = {}

        for plugin_name, plugin_record in registry.plugins.items():
            plugin_data = plugin_record.as_mapping()
            if plugin_data.get("category", "calculation").lower() != "calculation":
                # Ignore non-property plugins
                continue

            if self._model.loaded_from_process:
                is_installed = is_plugin_installed(plugin_name)
            else:
                status = (
                    statuses[plugin_name]
                    if statuses is not None
                    else get_plugin_status(
                        plugin_name,
                        plugin_record,
                        check_runtime=False,
                    )
                )
                is_installed = status.is_installed
                if status.is_incompatible:
                    self.incompatible_plugin_data[status.distribution_name] = {
                        **plugin_data,
                        "installed_version": status.installed_version,
                        "plugin_compatible": status.plugin_compatible,
                        "app_compatible": status.app_compatible,
                        "requirement_error": status.requirement_error,
                        "requirement": status.requirement,
                        "registration_error": status.registration_error,
                    }

            if not is_installed:
                self.available_plugins_list.append(plugin_data["title"])

        self._model.incompatible_plugins = {
            identifier
            for identifier, owner in self.entry_point_distributions.items()
            if owner and canonicalize_name(owner) in self.incompatible_plugin_data
        }

        self._model.available_plugins_fetched = True

    def _schedule_runtime_plugin_checks(self) -> None:
        if self._model.loaded_from_process:
            return

        loading_message = """
            <div class="loading" style="display: flex; align-items: center; font-size: unset;">
                Loading/validating plugins...
                <i class="fa fa-spinner fa-spin fa-1x fa-fw"></i>
            </div>
        """

        self.available_plugins.value = loading_message
        self.installed_plugins.children = [ipw.HTML(value=loading_message)]

        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            asyncio.run(self._refresh_runtime_plugin_checks())
        else:
            self._plugin_status_task = loop.create_task(
                self._refresh_runtime_plugin_checks()
            )

    async def _refresh_runtime_plugin_checks(self) -> None:
        registry = self._plugin_registry
        plugin_records = {
            name: record
            for name, record in registry.plugins.items()
            if record.category.lower() == "calculation"
        }
        statuses = await get_plugin_statuses(plugin_records)
        self._fetch_available_plugins(statuses=statuses)

        for identifier, model in self._model.get_models():
            if identifier in self._model.incompatible_plugins:
                model.include = False

        self.incompatible_plugins_list = [
            plugin_data["title"]
            for plugin_data in self.incompatible_plugin_data.values()
        ]
        self.installed_plugins_list = [
            row
            for identifier, row in self.installed_plugin_rows.items()
            if identifier not in self._model.incompatible_plugins
        ]
        self.installed_plugins.children = self.installed_plugins_list
        self._set_incompatible_plugins()
        self._toggle_incompatible_plugins()
        self._set_available_plugins()
