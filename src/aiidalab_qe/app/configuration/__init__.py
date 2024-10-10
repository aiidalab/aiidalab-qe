"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.panel import Panel, SettingPanel
from aiidalab_widgets_base import WizardAppWidgetStep

from .advanced import AdvancedSettings
from .model import ConfigurationModel
from .workflow import WorkChainSettings

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class ConfigureQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    def __init__(self, model: ConfigurationModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workflow configuration panel")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_confirmation_change,
            "confirmed",
        )
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.workchain.observe(
            self._on_protocol_change,
            "protocol",
        )

        self.missing_structure_message = """
            <div class="alert alert-info">
                <b>Please set the input structure first.</b>
            </div>
        """
        self.structure_set_message = ipw.HTML(self.missing_structure_message)

        self.workchain_settings = WorkChainSettings(config_model=model)
        self.advanced_settings = AdvancedSettings(config_model=model)

        self.built_in_settings = [
            self.workchain_settings,
            self.advanced_settings,
        ]

        self.settings: dict[str, Panel] = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        self._fetch_setting_entries()

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.tab = ipw.Tab(
            children=self.built_in_settings,
            layout=ipw.Layout(min_height="250px"),
            selected_index=None,
        )
        for i, setting in enumerate(self.built_in_settings):
            self.tab.set_title(i, setting.title)
        self.tab.observe(
            self._on_tab_change,
            "selected_index",
        )
        self.tab.selected_index = 0

        self.confirm_button = ipw.Button(
            description="Confirm",
            tooltip="Confirm the currently selected settings and go to the next step.",
            button_style="success",
            icon="check-circle",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        ipw.dlink(
            (self, "state"),
            (self.confirm_button, "disabled"),
            lambda state: state != self.State.CONFIGURED,
        )
        self.confirm_button.on_click(self.confirm)

        self.children = [
            self.structure_set_message,
            *self.property_children,
            self.tab,
            self.confirm_button,
        ]

        self.rendered = True

    def get_configuration_parameters(self):
        return self._model.get_model_state()

    def set_configuration_parameters(self, parameters):
        self._model.set_model_state(parameters)

    def is_saved(self):
        return self._model.confirmed

    def confirm(self, _=None):
        self._model.configuration_parameters = self.get_configuration_parameters()
        self._model.confirmed = True

    def reset(self):
        self._model.reset()
        if self.rendered:
            self.tab.selected_index = 0
        for _, settings in self.settings.items():
            settings.reset()

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update_state()

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: SettingPanel = self.tab.children[tab_index]  # type: ignore
        tab.render()
        tab.update()

    def _on_input_structure_change(self, _):
        self._update_missing_structure_warning()
        self.reset()

    def _on_protocol_change(self, _):
        self._model.advanced.update()

    def _on_confirmation_change(self, _):
        if not self._model.confirmed:
            self._model.configuration_parameters = {}
        self._update_state()

    def _fetch_setting_entries(self):
        """Handle plugin specific settings."""

        self.properties = {}
        self.reminder_info = {}
        self.property_children = [ipw.HTML("Select which properties to calculate:")]

        outlines = get_entry_items("aiidalab_qe.properties", "outline")
        models = get_entry_items("aiidalab_qe.properties", "model")
        settings = get_entry_items("aiidalab_qe.properties", "setting")
        for identifier in settings:
            model = models[identifier]()
            self._model.add_model(identifier, model)

            outline = outlines[identifier]()
            info = ipw.HTML()
            ipw.link(
                (model, "include"),
                (outline.include, "value"),
            )

            def toggle_plugin(change, identifier=identifier, info=info):
                if change["new"]:
                    info.value = f"Customize {identifier} settings below"
                else:
                    info.value = ""
                self._update_panel()

            model.observe(
                toggle_plugin,
                "include",
            )

            self.properties[identifier] = outline
            self.property_children.append(
                ipw.HBox(
                    children=[
                        outline,
                        info,
                    ]
                )
            )

            self.settings[identifier] = settings[identifier](
                parent=self,
                identifier=identifier,
                config_model=self._model,
            )

    def _update_missing_structure_warning(self):
        self.structure_set_message.value = (
            self.missing_structure_message
            if self._model.input_structure is None
            else ""
        )

    def _update_panel(self, _=None):
        self.tab.children = self.built_in_settings
        for identifier in self.properties:
            model = self._model.get_model(identifier)
            setting = self.settings[identifier]
            if model and model.include:
                self.tab.children += (setting,)
                self.tab.set_title(
                    len(self.tab.children) - 1,
                    setting.title,
                )

    def _update_state(self, _=None):
        if self._model.confirmed:
            self.state = self.State.SUCCESS
        elif self.previous_step_state is self.State.SUCCESS:
            self.state = self.State.CONFIGURED
        elif self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
        else:
            self.state = self.State.INIT
