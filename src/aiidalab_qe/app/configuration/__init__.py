"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_widgets_base import WizardAppWidgetStep

from .advanced import AdvancedSettings
from .workflow import WorkChainSettings


class ConfigureQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    confirmed = tl.Bool()
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)
    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    # output dictionary
    configuration_parameters = tl.Dict()

    def __init__(self, **kwargs):
        self.workchain_settings = WorkChainSettings()
        self.advanced_settings = AdvancedSettings()

        ipw.dlink(
            (self.workchain_settings.workchain_protocol, "value"),
            (self.advanced_settings, "protocol"),
        )
        ipw.dlink(
            (self.workchain_settings.spin_type, "value"),
            (self.advanced_settings, "spin_type"),
        )
        ipw.dlink(
            (self.workchain_settings.electronic_type, "value"),
            (self.advanced_settings, "electronic_type"),
        )
        ipw.dlink(
            (self, "input_structure"),
            (self.advanced_settings, "input_structure"),
        )
        #
        ipw.dlink(
            (self, "input_structure"),
            (self.workchain_settings, "input_structure"),
        )
        #
        self.built_in_settings = [
            self.workchain_settings,
            self.advanced_settings,
        ]
        self.tab = ipw.Tab(
            children=self.built_in_settings,
            layout=ipw.Layout(min_height="250px"),
        )

        self.tab.set_title(0, "Basic settings")
        self.tab.set_title(1, "Advanced settings")

        # store the property identifier and setting panel for all plugins
        # only show the setting panel when the corresponding property is selected
        # first add the built-in settings
        self.settings = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        # list of trailets to link
        # if new trailets are added to the settings, they need to be added here
        trailets_list = ["input_structure", "protocol", "electronic_type", "spin_type"]

        # then add plugin specific settings
        entries = get_entry_items("aiidalab_qe.properties", "setting")
        for identifier, entry_point in entries.items():
            self.settings[identifier] = entry_point(parent=self)
            self.settings[identifier].identifier = identifier
            # link basic protocol to all plugin specific protocols
            if identifier in self.workchain_settings.properties:
                self.workchain_settings.properties[identifier].run.observe(
                    self._update_panel, "value"
                )
            # link the trailets if they exist in the plugin specific settings
            for trailet in trailets_list:
                if hasattr(self.settings[identifier], trailet):
                    ipw.dlink(
                        (self.advanced_settings, trailet),
                        (self.settings[identifier], trailet),
                    )

        self._submission_blocker_messages = ipw.HTML()

        self.confirm_button = ipw.Button(
            description="Confirm",
            tooltip="Confirm the currently selected settings and go to the next step.",
            button_style="success",
            icon="check-circle",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )

        self.confirm_button.on_click(self.confirm)

        super().__init__(
            children=[
                self.tab,
                self._submission_blocker_messages,
                self.confirm_button,
            ],
            **kwargs,
        )

    @tl.observe("previous_step_state")
    def _observe_previous_step_state(self, _change):
        self._update_state()

    @tl.observe("input_structure")
    def _observe_input_structure(self, _change):
        if self.input_structure is not None and self.input_structure.pbc == (False, False, False): 
            pass

    def get_configuration_parameters(self):
        """Get the parameters of the configuration step."""

        return {s.identifier: s.get_panel_value() for s in self.tab.children}

    def set_configuration_parameters(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        with self.hold_trait_notifications():
            for identifier, settings in self.settings.items():
                if parameters.get(identifier):
                    settings.set_panel_value(parameters[identifier])

    def _update_state(self, _=None):
        if self.previous_step_state == self.State.SUCCESS:
            self.confirm_button.disabled = False
            self._submission_blocker_messages.value = ""
            self.state = self.State.CONFIGURED
            # update plugin specific settings
            for _, settings in self.settings.items():
                settings._update_state()
        elif self.previous_step_state == self.State.FAIL:
            self.state = self.State.FAIL
        else:
            self.confirm_button.disabled = True
            self.state = self.State.INIT
            self.reset()

    def confirm(self, _=None):
        self.configuration_parameters = self.get_configuration_parameters()
        self.confirm_button.disabled = False
        self.state = self.State.SUCCESS

    def is_saved(self):
        """Check if the current step is saved.
        That all changes are confirmed.
        """
        new_parameters = self.get_configuration_parameters()
        return new_parameters == self.configuration_parameters

    @tl.default("state")
    def _default_state(self):
        return self.State.INIT
    
    def reset(self):
        """Reset the widgets in all settings to their initial states."""
        with self.hold_trait_notifications():
            self.input_structure = None
            for _, settings in self.settings.items():
                settings.reset()

    def _update_panel(self, _=None):
        """Dynamic add/remove the panel based on the selected properties."""
        # only keep basic and advanced settings
        self.tab.children = self.built_in_settings
        # add plugin specific settings
        for identifier in self.workchain_settings.properties:
            if (
                identifier in self.settings
                and self.workchain_settings.properties[identifier].run.value
            ):
                self.tab.children += (self.settings[identifier],)
                self.tab.set_title(
                    len(self.tab.children) - 1, self.settings[identifier].title
                )
