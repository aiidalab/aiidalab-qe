# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""
from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiidalab_widgets_base import WizardAppWidgetStep

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items

from .advanced import AdvancedSettings
from .workflow import WorkChainSettings


class ConfigureQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    confirmed = tl.Bool()
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)
    workchain_settings = tl.Instance(WorkChainSettings, allow_none=True)
    advanced_settings = tl.Instance(AdvancedSettings, allow_none=True)
    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    # output dictionary
    configuration_parameters = tl.Dict()

    def __init__(self, **kwargs):
        self.workchain_settings = WorkChainSettings()
        self.workchain_settings.relax_type.observe(self._update_state, "value")
        self.advanced_settings = AdvancedSettings()

        ipw.dlink(
            (self.workchain_settings.workchain_protocol, "value"),
            (self.advanced_settings, "protocol"),
        )
        ipw.dlink(
            (self, "input_structure"),
            (self.advanced_settings, "input_structure"),
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

        self.tab.set_title(0, "Workflow")
        self.tab.set_title(1, "Advanced settings")

        # add plugin specific settings
        self.settings = {}
        # add plugin specific settings
        entries = get_entry_items("aiidalab_qe.properties", "setting")
        for name, entry_point in entries.items():
            self.settings[name] = entry_point(parent=self)
            self.settings[name].identifier = name
            # link basic protocol to all plugin specific protocols
            if hasattr(self.settings[name], "workchain_protocol"):
                ipw.dlink(
                    (self.workchain_settings.workchain_protocol, "value"),
                    (self.settings[name].workchain_protocol, "value"),
                )
            if name in self.workchain_settings.properties:
                self.workchain_settings.properties[name].run.observe(
                    self._update_panel, "value"
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
    def _observe_previous_step_state(self, change):
        self._update_state()

    def get_configuration_parameters(self):
        """Get the parameters of the configuration step."""
        # create the the initial_magnetic_moments as None (Default)
        initial_magnetic_moments = None
        # create the override parameters for sub PwBaseWorkChain
        pw_overrides = {"base": {}, "scf": {}, "nscf": {}, "band": {}}
        for key in ["base", "scf", "nscf", "band"]:
            pw_overrides[key][
                "pseudo_family"
            ] = self.advanced_settings.pseudo_family_selector.value

            if self.advanced_settings.pseudo_setter.pseudos:
                pw_overrides[key].setdefault("pw", {"parameters": {"SYSTEM": {}}})
                pw_overrides[key]["pw"][
                    "pseudos"
                ] = self.advanced_settings.pseudo_setter.pseudos

                pw_overrides[key]["pw"]["parameters"]["SYSTEM"][
                    "ecutwfc"
                ] = self.advanced_settings.pseudo_setter.ecutwfc
                pw_overrides[key]["pw"]["parameters"]["SYSTEM"][
                    "ecutrho"
                ] = self.advanced_settings.pseudo_setter.ecutrho

            if self.advanced_settings.override.value:
                pw_overrides[key].setdefault("pw", {"parameters": {"SYSTEM": {}}})
                pw_overrides[key]["pw"]["parameters"]["SYSTEM"][
                    "tot_charge"
                ] = self.advanced_settings.value.get("total_charge")
                if self.workchain_settings.spin_type.value == "collinear":
                    initial_magnetic_moments = (
                        self.advanced_settings.magnetization.get_magnetization()
                    )

                if key in ["base", "scf"]:
                    pw_overrides[key][
                        "kpoints_distance"
                    ] = self.advanced_settings.value.get("kpoints_distance")

                    if self.workchain_settings.electronic_type.value == "metal":
                        # smearing type setting
                        pw_overrides[key]["pw"]["parameters"]["SYSTEM"][
                            "smearing"
                        ] = self.advanced_settings.value.get("smearing")

                        # smearing degauss setting
                        pw_overrides[key]["pw"]["parameters"]["SYSTEM"][
                            "degauss"
                        ] = self.advanced_settings.value.get("degauss")

        overrides = {
            "relax": {
                "base": pw_overrides["base"],
            },
            "bands": {
                "scf": pw_overrides["scf"],
                "bands": pw_overrides["band"],
            },
            "pdos": {
                "scf": pw_overrides["scf"],
                "nscf": pw_overrides["nscf"],
            },
        }
        # Work chain settings
        workchain_settings = self.workchain_settings.get_setting_parameters()
        return {
            "workchain_settings": workchain_settings,
            "overrides": overrides,
            "initial_magnetic_moments": initial_magnetic_moments,
        }

    def set_input_parameters(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        with self.hold_trait_notifications():
            # Work chain settings
            workchain_settings = parameters.get("workchain_settings", {})
            self.workchain_settings.set_setting_parameters(workchain_settings)
            # Advanced settings
            if parameters.get("pseudo_family", False):
                self.advanced_settings.pseudo_family_selector.value = parameters[
                    "pseudo_family"
                ]
            if parameters.get("kpoints_distance_override", None) is not None:
                self.advanced_settings.kpoints.distance.value = parameters[
                    "kpoints_distance_override"
                ]
                self.advanced_settings.kpoints.override.value = True
            if parameters.get("degauss_override", None) is not None:
                self.advanced_settings.smearing.degauss.value = parameters[
                    "degauss_override"
                ]
                self.advanced_settings.smearing.override.value = True
            if parameters.get("smearing_override", None) is not None:
                self.advanced_settings.smearing.smearing.value = parameters[
                    "smearing_override"
                ]
                self.advanced_settings.smearing.override.value = True

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
            self.set_input_parameters(DEFAULT_PARAMETERS)

    def confirm(self, _=None):
        self.configuration_parameters = self.get_configuration_parameters()
        self.confirm_button.disabled = False
        self.state = self.State.SUCCESS

    @tl.default("state")
    def _default_state(self):
        return self.State.INIT

    def reset(self):
        with self.hold_trait_notifications():
            self.set_input_parameters(DEFAULT_PARAMETERS)

    def _update_panel(self, _=None):
        """Dynamic add/remove the panel based on the selected properties."""
        # only keep basic and advanced settings
        self.tab.children = self.built_in_settings
        # add plugin specific settings
        for name in self.workchain_settings.properties:
            if (
                name in self.settings
                and self.workchain_settings.properties[name].run.value
            ):
                self.tab.children += (self.settings[name],)
                self.tab.set_title(
                    len(self.tab.children) - 1, self.settings[name].title
                )
