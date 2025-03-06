"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import (
    ConfigurationSettingsModel,
    ConfigurationSettingsPanel,
)
from aiidalab_qe.common.widgets import LinkButton
from aiidalab_qe.common.wizard import QeConfirmableDependentWizardStep

from .advanced import (
    AdvancedConfigurationSettingsModel,
    AdvancedConfigurationSettingsPanel,
)
from .basic import BasicConfigurationSettingsModel, BasicConfigurationSettingsPanel
from .model import ConfigurationStepModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class ConfigureQeAppWorkChainStep(
    QeConfirmableDependentWizardStep[ConfigurationStepModel]
):
    missing_information_warning = "Missing input structure. Please set it first."

    def __init__(self, model: ConfigurationStepModel, **kwargs):
        super().__init__(
            model=model,
            confirm_kwargs={
                "tooltip": "Confirm the currently selected settings and go to the next step",
            },
            **kwargs,
        )
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
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
        # orbit is set from loaded process
        ipw.dlink(
            (self._model.get_model("advanced"), "spin_orbit"),
            (self._model.get_model("workchain"), "spin_orbit"),
        )

        self.settings = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        self.installed_property_children = []
        self.not_installed_property_children = []
        self._fetch_plugin_calculation_settings()

    def _render(self):
        super()._render()

        self._fetch_not_installed_property()

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
            button_style="primary",  # Keeps it prominent
            tooltip="Browse and install additional plugins from the Plugin Store",
        )

        self.sub_steps = ipw.Accordion(
            children=[
                ipw.VBox(
                    children=[
                        InAppGuide(identifier="properties-selection"),
                        *self.installed_property_children,
                        ipw.HTML("<hr>"),
                        ipw.HTML(
                            value="""
                            <p style="font-size:14px; line-height:1.6;">
                                The following additional property calculations are available in the
                                <b>Plugin registry</b> but are currently disabled because the required plugins are not installed.
                            </p>
                            <p style="font-size:14px; line-height:1.6;">
                                To enable them, please visit the <b>Plugin store</b> and install the necessary plugins.
                            </p>
                            """,
                            layout=ipw.Layout(margin="10px 0px"),
                        ),
                        self.install_new_plugin_button,
                        *self.not_installed_property_children,
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
        ]

        self.children = [
            self.content,
            self.confirm_box,
        ]

    def _post_render(self):
        self._update_tabs()

    def reset(self):
        self._model.reset()
        if self.rendered:
            self.sub_steps.selected_index = None
            self.tabs.selected_index = 0

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: ConfigurationSettingsPanel = self.tabs.children[tab_index]  # type: ignore
        tab.render()
        tab.update()

    def _on_input_structure_change(self, _):
        self._model.update()
        self.reset()

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

    def _update_state(self, _=None):
        if self._model.confirmed:
            self.state = self.State.SUCCESS
        elif self.previous_step_state is self.State.SUCCESS:
            self.state = self.State.CONFIGURED
        elif self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
        else:
            self.state = self.State.INIT

    def _fetch_not_installed_property(self, plugin_config_source=None):
        from aiidalab_qe.app.utils.plugin_manager import (
            DEFAULT_PLUGIN_CONFIG_SOURCE,
            PluginManager,
            is_package_installed,
        )

        self.not_installed_property_children = []

        plugin_config_source = plugin_config_source or DEFAULT_PLUGIN_CONFIG_SOURCE
        plugin_manager = PluginManager(plugin_config_source)
        for plugin_name, plugin_data in plugin_manager.data.items():
            if (
                plugin_data.get("category", "calculation").lower() != "calculation"
            ):  # Ignore non-property plugins
                continue

            is_installed = is_package_installed(plugin_name)
            if not is_installed:
                checkbox = ipw.Checkbox(
                    value=False,
                    description=plugin_data["title"],
                    disabled=True,
                    indent=False,
                    style={"description_width": "initial"},
                )
                self.not_installed_property_children.append(checkbox)

    def _fetch_plugin_calculation_settings(self):
        outlines = get_entry_items("aiidalab_qe.properties", "outline")
        entries = get_entry_items("aiidalab_qe.properties", "configuration")
        for identifier, configuration in entries.items():
            for key in ("panel", "model"):
                if key not in configuration:
                    raise ValueError(f"Entry {identifier} is missing the '{key}' key")

            model: ConfigurationSettingsModel = configuration["model"]()
            self._model.add_model(identifier, model)

            outline = outlines[identifier]()
            info = ipw.HTML()
            ipw.link(
                (model, "include"),
                (outline.include, "value"),
            )

            if identifier == "bands":
                ipw.dlink(
                    (self._model, "input_structure"),
                    (outline.include, "disabled"),
                    lambda _: not self._model.has_pbc,
                )

            def toggle_plugin(change, identifier=identifier, model=model, info=info):
                if change["new"]:
                    info.value = (
                        f"Customize {identifier} settings in <b>Step 2.2</b> if needed"
                    )
                else:
                    info.value = ""
                model.update()
                self._update_tabs()

            model.observe(
                toggle_plugin,
                "include",
            )

            self.installed_property_children.append(
                ipw.HBox(
                    children=[
                        outline,
                        info,
                    ]
                )
            )

            panel: ConfigurationSettingsPanel = configuration["panel"](model=model)
            self.settings[identifier] = panel
