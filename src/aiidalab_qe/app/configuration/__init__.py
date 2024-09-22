"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiida_quantumespresso.common.types import RelaxType
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_widgets_base import WizardAppWidgetStep

from .advanced import AdvancedSettings
from .model import AdvancedModel, ConfigurationModel, WorkChainModel
from .workflow import WorkChainSettings

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class ConfigureQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)
    confirmed = tl.Bool()

    _no_structure_warning = """
        <div style="color: red;">
            Please set the input structure first.
        </div>
    """

    def __init__(self, model: ConfigurationModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workflow configuration panel")],
            **kwargs,
        )

        self._model = model

        workchain_model = WorkChainModel()
        advanced_model = AdvancedModel()
        ipw.dlink(
            (self._model, "input_structure"),
            (advanced_model, "input_structure"),
        )
        ipw.dlink(
            (workchain_model, "protocol"),
            (advanced_model, "protocol"),
        )
        ipw.dlink(
            (workchain_model, "spin_type"),
            (advanced_model, "spin_type"),
        )
        ipw.dlink(
            (workchain_model, "electronic_type"),
            (advanced_model, "electronic_type"),
        )
        # TODO necessary?
        self._model.observe(
            lambda _: advanced_model.update_kpoints_mesh(),
            "input_structure",
        )

        self.workchain_settings = WorkChainSettings(model=workchain_model)
        self.advanced_settings = AdvancedSettings(model=advanced_model)

        self.built_in_settings = [
            self.workchain_settings,
            self.advanced_settings,
        ]

        self.settings = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        self._fetch_setting_entries()

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.structure_set_message = ipw.HTML()
        ipw.dlink(
            (self._model, "input_structure"),
            (self.structure_set_message, "value"),
            lambda structure: self._no_structure_warning if structure is None else "",
        )

        self.tab = ipw.Tab(
            children=self.built_in_settings,
            layout=ipw.Layout(min_height="250px"),
            selected_index=None,
        )
        self.tab.set_title(0, "Basic settings")
        self.tab.set_title(1, "Advanced settings")
        self.tab.observe(self._on_tab_change, "selected_index")
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

    def reset(self):
        """Reset the widgets in all settings to their initial states."""
        with self.hold_trait_notifications():
            for _, settings in self.settings.items():
                if settings.rendered:
                    settings.reset()
            for key, p in self.properties.items():
                p.run.value = key in DEFAULT["workchain"]["properties"]

    def is_saved(self):
        """Check if the current step is saved.
        That all changes are confirmed.
        """
        new_parameters = self.get_configuration_parameters()
        return new_parameters == self._model.configuration_parameters

    def confirm(self, _=None):
        self._model.configuration_parameters = self.get_configuration_parameters()
        self.confirm_button.disabled = False
        self.state = self.State.SUCCESS

    def get_configuration_parameters(self):
        parameters = {
            setting.identifier: setting.get_panel_value()
            for setting in self.settings.values()
        }
        properties = self._get_properties()
        # TODO necessary to store the properties in the workchain settings?
        parameters["workchain"].update("properties", properties)

    def set_configuration_parameters(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        with self.hold_trait_notifications():
            for identifier, settings in self.settings.items():
                if parameters.get(identifier):
                    settings.set_panel_value(parameters[identifier])
            properties = parameters.get("properties", [])
            for name in self.properties:
                if name in properties:
                    self.properties[name].run.value = True
                else:
                    self.properties[name].run.value = False

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, change):
        self._update_state(change["new"])

    def _on_tab_change(self, change):
        if (tab := change["new"]) is None:
            return
        self.tab.children[tab].render()  # type: ignore

    def _fetch_setting_entries(self):
        """Handle plugin specific settings."""

        self.properties = {}
        self.reminder_info = {}
        self.property_children = [ipw.HTML("Select which properties to calculate:")]
        entries = get_entry_items("aiidalab_qe.properties", "outline")
        models = get_entry_items("aiidalab_qe.properties", "model")
        settings = get_entry_items("aiidalab_qe.properties", "setting")
        for (name, entry_point), setting in zip(entries.items(), settings.values()):
            self.properties[name] = entry_point()
            self.properties[name].run.observe(self._update_panel, "value")
            self.reminder_info[name] = ipw.HTML()
            self.property_children.append(
                ipw.HBox(
                    children=[
                        self.properties[name],
                        self.reminder_info[name],
                    ]
                )
            )

            def update_reminder_info(change, name=name):
                info = self.reminder_info[name]
                info.value = (
                    f"Customize {name} settings in the corresponding tab"
                    if change["new"]
                    else ""
                )

            if name in settings:
                self.properties[name].run.observe(update_reminder_info, "value")
                kwargs = {"parent": self, "identifier": name}
                # TODO drop check in the future - plugin models should be required
                if name in models:
                    kwargs["model"] = models[name]()
                self.settings[name] = setting(**kwargs)

        self.property_children.append(
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
                    The band structure workflow will automatically detect the default
                    path in reciprocal space using the
                    <a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">SeeK-path tool</a>.
                </div>
            """)
        )

    def _update_panel(self, _=None):
        """Dynamic add/remove the panel based on the selected properties."""
        self.tab.children = self.built_in_settings
        for identifier in self.properties:
            if identifier in self.settings and self.properties[identifier].run.value:
                self.tab.children += (self.settings[identifier],)
                self.tab.set_title(
                    len(self.tab.children) - 1, self.settings[identifier].title
                )

    def _update_state(self, previous_step_state):
        if previous_step_state == self.State.SUCCESS:
            self.state = self.State.CONFIGURED
            for settings in self.settings.values():
                settings._update_state()
        elif previous_step_state == self.State.FAIL:
            self.state = self.State.FAIL
        else:
            self.state = self.State.INIT
            self.reset()

    def _get_properties(self):
        properties = []
        run_bands = False
        run_pdos = False
        for name in self.properties:
            if self.properties[name].run.value:
                properties.append(name)
            if name == "bands":
                run_bands = True
            elif name == "pdos":
                run_bands = True

        if RelaxType(self._model.basic.relax_type) is not RelaxType.NONE or not (
            run_bands or run_pdos
        ):
            properties.append("relax")
        return properties
