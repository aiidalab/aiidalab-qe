"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.app.utils import get_entry_items
from aiidalab_widgets_base import WizardAppWidgetStep

from .advanced import AdvancedSettings
from .model import config_model as model
from .workflow import WorkChainSettings


class ConfigureQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)
    confirmed = tl.Bool()

    _structure_not_set_warning = """
        <div style="color: red;">
            Please set the input structure first.
        </div>
    """

    def __init__(self, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workflow configuration panel")],
            **kwargs,
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.structure_set_message = ipw.HTML()
        ipw.dlink(
            (model, "input_structure"),
            (self.structure_set_message, "value"),
            lambda s: "" if s else self._structure_not_set_warning,
        )

        self.workchain_settings = WorkChainSettings(callback=self._update_panel)
        self.advanced_settings = AdvancedSettings()

        self.built_in_settings = [
            self.workchain_settings,
            self.advanced_settings,
        ]
        self.tab = ipw.Tab(
            children=self.built_in_settings,
            layout=ipw.Layout(min_height="250px"),
            selected_index=None,
        )
        self.tab.set_title(0, "Basic settings")
        self.tab.set_title(1, "Advanced settings")

        self.tab.observe(self._on_tab_change, "selected_index")

        self.tab.selected_index = 0

        # store the property identifier and setting panel for all plugins
        # only show the setting panel when the corresponding property is selected
        # first add the built-in settings
        self.settings = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        # then add plugin specific settings
        entries = get_entry_items("aiidalab_qe.properties", "setting")
        for identifier, entry_point in entries.items():
            self.settings[identifier] = entry_point(parent=self)
            self.settings[identifier].identifier = identifier

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
            self.tab,
            self.confirm_button,
        ]

        self.rendered = True

    def reset(self):
        """Reset the widgets in all settings to their initial states."""
        with self.hold_trait_notifications():
            for _, settings in self.settings.items():
                if settings.rendered:
                    settings.reset()

    def is_saved(self):
        """Check if the current step is saved.
        That all changes are confirmed.
        """
        new_parameters = self.get_configuration_parameters()
        return new_parameters == model.configuration_parameters

    def confirm(self, _=None):
        model.configuration_parameters = self.get_configuration_parameters()
        self.confirm_button.disabled = False
        self.state = self.State.SUCCESS

    def get_configuration_parameters(self):
        """Get the parameters of the configuration step."""
        if not hasattr(self, "tab"):
            return {}  # TODO restructure plugin system to use models
        return {s.identifier: s.get_panel_value() for s in self.tab.children}

    def set_configuration_parameters(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        with self.hold_trait_notifications():
            for identifier, settings in self.settings.items():
                if parameters.get(identifier):
                    settings.set_panel_value(parameters[identifier])

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update_state()

    def _update_state(self):
        if self.previous_step_state == self.State.SUCCESS:
            self.state = self.State.CONFIGURED
            # for settings in self.settings.values():
            #     settings._update_state()
        elif self.previous_step_state == self.State.FAIL:
            self.state = self.State.FAIL
        else:
            self.state = self.State.INIT
            # self.reset()  # TODO redundant?

    def _on_tab_change(self, change):
        if (tab := change["new"]) is None:
            return
        self.tab.children[tab].render()  # type: ignore

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
