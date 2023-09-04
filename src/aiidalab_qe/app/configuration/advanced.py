# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""
import os

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from IPython.display import clear_output, display

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS


class AdvancedSettings(ipw.VBox):
    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Advanced Settings</h4></div>"""
    )
    description = ipw.HTML("""Select the advanced settings for the <b>pw.x</b> code.""")
    kpoints_description = ipw.HTML(
        """<div>
        The k-points mesh density of the SCF calculation is set by the <b>protocol</b>.
        The value below represents the maximum distance between the k-points in each direction of reciprocal space.
        Tick the box to override the default, smaller is more accurate and costly. </div>"""
    )

    # protocol interface
    protocol = tl.Unicode(allow_none=True)

    # output dictionary
    value = tl.Dict()

    def __init__(self, default_protocol=None, **kwargs):
        self._default_protocol = default_protocol or DEFAULT_PARAMETERS["protocol"]

        # Override setting widget
        self.override_prompt = ipw.HTML("<b>&nbsp;&nbsp;Override&nbsp;</b>")
        self.override = ipw.Checkbox(
            description="",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="10%"),
        )
        self.override_widget = ipw.HBox(
            [self.override_prompt, self.override],
            layout=ipw.Layout(max_width="20%"),
        )
        # Smearing setting widget
        self.smearing = SmearingSettings()
        ipw.dlink(
            (self.override, "value"),
            (self.smearing, "disabled"),
            lambda override: not override,
        )
        self.smearing.observe(
            self._callback_value_set, ["degauss_value", "smearing_value"]
        )

        # Kpoints setting widget
        self.kpoints_distance = ipw.FloatText(
            step=0.05,
            description="K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override, "value"),
            (self.kpoints_distance, "disabled"),
            lambda override: not override,
        )
        self.kpoints_distance.observe(self._callback_value_set, "value")

        # Total change setting widget
        self.total_charge = ipw.BoundedFloatText(
            min=-3,
            max=3,
            step=0.01,
            disabled=False,
            description="Total charge:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override, "value"),
            (self.total_charge, "disabled"),
            lambda override: not override,
        )
        self.total_charge.observe(self._callback_value_set, "value")

        self.magnetization = MagnetizationSettings()
        self.list_overrides = [
            self.magnetization.override,
        ]
        for override in self.list_overrides:
            ipw.dlink(
                (self.override, "value"),
                (override, "disabled"),
                lambda override: not override,
            )
        self.override.observe(self.set_advanced_settings, "value")
        super().__init__(
            children=[
                self.title,
                ipw.HBox(
                    [self.description, self.override_widget],
                    layout=ipw.Layout(height="50px", justify_content="space-between"),
                ),
                # total charge setting widget
                self.total_charge,
                # magnetization setting widget
                self.magnetization,
                # smearing setting widget
                self.smearing,
                # Kpoints setting widget
                self.kpoints_description,
                self.kpoints_distance,
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

        # Default settings to trigger the callback
        self.reset()

    @tl.observe("protocol")
    def _protocol_changed(self, _):
        """Input protocol changed, update the widget values."""
        self._update_settings_from_protocol(self.protocol)

    def _update_settings_from_protocol(self, protocol):
        """Update the values of sub-widgets from the given protocol, this will
        trigger the callback of the sub-widget if it is exist.
        """
        self.smearing.protocol = protocol

        parameters = PwBaseWorkChain.get_protocol_inputs(protocol)

        self.kpoints_distance.value = parameters["kpoints_distance"]

    def _callback_value_set(self, _=None):
        """Callback function to set the parameters"""
        settings = {
            "kpoints_distance": self.kpoints_distance.value,
            "total_charge": self.total_charge.value,
            "degauss": self.smearing.degauss_value,
            "smearing": self.smearing.smearing_value,
        }

        self.update_settings(**settings)

    def update_settings(self, **kwargs):
        """Set the output dict from the given keyword arguments.
        This function will only update the traitlets but not the widget value.

        This function can also be used to set values directly for testing purpose.
        """
        self.value = {k: v for k, v in kwargs.items()}

    def set_advanced_settings(self, _=None):
        self.magnetization.reset()

    def reset(self):
        """Reset the widget and the traitlets"""
        self.protocol = self._default_protocol

        with self.hold_trait_notifications():
            # Reset protocol dependent settings
            self._update_settings_from_protocol(self.protocol)
            # reset total charge
            self.total_charge.value = DEFAULT_PARAMETERS["tot_charge"]

            # reset the override checkbox
            self.override.value = False


class MagnetizationSettings(ipw.VBox):
    """Widget to set the initial magnetic moments for each kind names defined in the StructureData (StructureDtaa.get_kind_names())
    Usually these are the names of the elements in the StructureData
    (For example 'C' , 'N' , 'Fe' . However the StructureData can have defined kinds like 'Fe1' and 'Fe2')

    The widget generate a dictionary that can be used to set initial_magnetic_moments in the builder of PwBaseWorkChain

    Attributes:
        input_structure(StructureData): trait that containes the input_strucgure (confirmed structure from previous step)
    """

    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    def __init__(self, **kwargs):
        self.input_structure = orm.StructureData()
        self.input_structure_labels = []
        self.description = ipw.HTML(
            "Define magnetization: Input structure not confirmed"
        )
        self.kinds = self.create_kinds_widget()
        self.kinds_widget_out = ipw.Output()
        self.override = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        super().__init__(
            children=[
                ipw.HBox(
                    [
                        self.override,
                        self.description,
                        self.kinds_widget_out,
                    ],
                ),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )
        self.display_kinds()
        self.override.observe(self._disable_kinds_widgets, "value")

    def _disable_kinds_widgets(self, _=None):
        for i in range(len(self.kinds.children)):
            self.kinds.children[i].disabled = not self.override.value

    def reset(self):
        self.override.value = False
        if hasattr(self.kinds, "children") and self.kinds.children:
            for i in range(len(self.kinds.children)):
                self.kinds.children[i].value = 0.0

    def create_kinds_widget(self):
        if self.input_structure_labels:
            widgets_list = []
            for kind_label in self.input_structure_labels:
                kind_widget = ipw.BoundedFloatText(
                    description=kind_label,
                    min=-4,
                    max=4,
                    step=0.1,
                    value=0.0,
                    disabled=True,
                )
                widgets_list.append(kind_widget)
            kinds_widget = ipw.VBox(widgets_list)
        else:
            kinds_widget = None

        return kinds_widget

    def update_kinds_widget(self):
        self.input_structure_labels = self.input_structure.get_kind_names()
        self.kinds = self.create_kinds_widget()
        self.description.value = "Define magnetization: "
        self.display_kinds()

    def display_kinds(self):
        if "PYTEST_CURRENT_TEST" not in os.environ and self.kinds:
            with self.kinds_widget_out:
                clear_output()
                display(self.kinds)

    def _update_widget(self, change):
        self.input_structure = change["new"]
        self.update_kinds_widget()

    def get_magnetization(self):
        """Method to generate the dictionary with the initial magnetic moments"""
        magnetization = {}
        for i in range(len(self.kinds.children)):
            magnetization[self.input_structure_labels[i]] = self.kinds.children[i].value
        return magnetization

    def _set_magnetization_values(self, **kwargs):
        """Update used for conftest setting all magnetization to a value"""
        self.override.value = True
        with self.hold_trait_notifications():
            if "initial_magnetic_moments" in kwargs:
                for i in range(len(self.kinds.children)):
                    self.kinds.children[i].value = kwargs["initial_magnetic_moments"]


class SmearingSettings(ipw.VBox):
    # accept protocol as input and set the values
    protocol = tl.Unicode(allow_none=True)

    # The output of the widget is a dictionary with the values of smearing and degauss
    degauss_value = tl.Float()
    smearing_value = tl.Unicode()

    smearing_description = ipw.HTML(
        """<p>
        The smearing type and width is set by the chosen <b>protocol</b>.
        Tick the box to override the default, not advised unless you've mastered <b>smearing effects</b> (click <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
        target="_blank">here</a> for a discussion).
    </p>"""
    )
    disabled = tl.Bool()

    def __init__(self, default_protocol=None, **kwargs):
        self._default_protocol = default_protocol or DEFAULT_PARAMETERS["protocol"]

        self.smearing = ipw.Dropdown(
            options=["cold", "gaussian", "fermi-dirac", "methfessel-paxton"],
            description="Smearing type:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.degauss = ipw.FloatText(
            step=0.005,
            description="Smearing width (Ry):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self, "disabled"),
            (self.degauss, "disabled"),
        )
        ipw.dlink(
            (self, "disabled"),
            (self.smearing, "disabled"),
        )
        self.degauss.observe(self._callback_value_set, "value")
        self.smearing.observe(self._callback_value_set, "value")

        super().__init__(
            children=[
                self.smearing_description,
                ipw.HBox([self.smearing, self.degauss]),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

        # Default settings to trigger the callback
        self.protocol = self._default_protocol

    @tl.default("disabled")
    def _default_disabled(self):
        return False

    @tl.observe("protocol")
    def _protocol_changed(self, _):
        """Input protocol changed, update the widget values."""
        self._update_settings_from_protocol(self.protocol)

    def _update_settings_from_protocol(self, protocol):
        """Update the widget values from the given protocol, and trigger the callback."""
        parameters = PwBaseWorkChain.get_protocol_inputs(protocol)["pw"]["parameters"][
            "SYSTEM"
        ]

        with self.hold_trait_notifications():
            # This changes will trigger callbacks
            self.degauss.value = parameters["degauss"]
            self.smearing.value = parameters["smearing"]

    def _callback_value_set(self, _=None):
        """callback function to set the smearing and degauss values"""
        settings = {
            "degauss": self.degauss.value,
            "smearing": self.smearing.value,
        }
        self.update_settings(**settings)

    def update_settings(self, **kwargs):
        """Set the output dict from the given keyword arguments.
        This function will only update the traitlets but not the widget value.
        """
        self.degauss_value = kwargs.get("degauss")
        self.smearing_value = kwargs.get("smearing")

    def reset(self):
        """Reset the widget and the traitlets"""
        self.protocol = self._default_protocol

        with self.hold_trait_notifications():
            self._update_settings_from_protocol(self.protocol)
            self.disabled = False
