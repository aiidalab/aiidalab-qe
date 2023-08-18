# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""
import os

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from IPython.display import clear_output, display


class AdvancedSettings(ipw.VBox):
    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Advanced Settings</h4></div>"""
    )
    description = ipw.HTML("""Select the advanced settings for the <b>pw.x</b> code.""")

    def __init__(self, **kwargs):
        self.override = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        self.smearing = SmearingSettings()
        self.kpoints = KpointSettings()
        self.tot_charge = TotalCharge()
        self.magnetization = MagnetizationSettings()
        self.list_overrides = [
            self.smearing.override,
            self.kpoints.override,
            self.tot_charge.override,
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
                    [
                        self.description,
                        self.override,
                    ],
                ),
                self.tot_charge,
                self.magnetization,
                self.smearing,
                self.kpoints,
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    def set_advanced_settings(self, _=None):
        self.smearing.reset()
        self.kpoints.reset()
        self.tot_charge.reset()
        self.magnetization.reset()


class TotalCharge(ipw.VBox):
    """Widget to define the total charge of the simulation"""

    tot_charge_default = tl.Float(default_value=0.0)

    def __init__(self, **kwargs):
        self.override = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        self.charge = ipw.BoundedFloatText(
            value=0,
            min=-3,
            max=3,
            step=0.01,
            disabled=False,
            description="Total charge:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override, "value"),
            (self.charge, "disabled"),
            lambda override: not override,
        )
        super().__init__(
            children=[
                ipw.HBox(
                    [
                        self.override,
                        self.charge,
                    ],
                ),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )
        self.charge.observe(self.set_tot_charge, "value")
        self.override.observe(self.set_tot_charge, "value")

    def set_tot_charge(self, _=None):
        self.charge.value = (
            self.charge.value if self.override.value else self.tot_charge_default
        )

    def _update_settings(self, **kwargs):
        """Update the override and override_tot_charge and override_tot_charge values by the given keyword arguments
        Therefore the override checkbox is not updated and defaults to True"""
        self.override.value = True
        with self.hold_trait_notifications():
            if "tot_charge" in kwargs:
                self.charge.value = kwargs["tot_charge"]

    def reset(self):
        with self.hold_trait_notifications():
            self.charge.value = self.tot_charge_default
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
                    min=-1,
                    max=1,
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
    smearing_description = ipw.HTML(
        """<p>
        The smearing type and width is set by the chosen <b>protocol</b>.
        Tick the box to override the default, not advised unless you've mastered <b>smearing effects</b> (click <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
        target="_blank">here</a> for a discussion).
    </p>"""
    )

    # The default of `smearing` and `degauss` the type and width
    # must be linked to the `protocol`
    degauss_default = tl.Float(default_value=0.01)
    smearing_default = tl.Unicode(default_value="cold")

    def __init__(self, **kwargs):
        self.override = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        self.smearing = ipw.Dropdown(
            options=["cold", "gaussian", "fermi-dirac", "methfessel-paxton"],
            value=self.smearing_default,
            description="Smearing type:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.degauss = ipw.FloatText(
            value=self.degauss_default,
            step=0.005,
            description="Smearing width (Ry):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override, "value"),
            (self.degauss, "disabled"),
            lambda override: not override,
        )
        ipw.dlink(
            (self.override, "value"),
            (self.smearing, "disabled"),
            lambda override: not override,
        )
        self.degauss.observe(self.set_smearing, "value")
        self.smearing.observe(self.set_smearing, "value")
        self.override.observe(self.set_smearing, "value")

        super().__init__(
            children=[
                self.smearing_description,
                ipw.HBox([self.override, self.smearing, self.degauss]),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    def set_smearing(self, _=None):
        self.degauss.value = (
            self.degauss.value if self.override.value else self.degauss_default
        )
        self.smearing.value = (
            self.smearing.value if self.override.value else self.smearing_default
        )

    def _update_settings(self, **kwargs):
        """Update the smearing and degauss values by the given keyword arguments
        This is the same as the `set_smearing` method but without the observer.
        Therefore the override checkbox is not updated and defaults to True"""
        self.override.value = True

        with self.hold_trait_notifications():
            if "smearing" in kwargs:
                self.smearing.value = kwargs["smearing"]

            if "degauss" in kwargs:
                self.degauss.value = kwargs["degauss"]

    def reset(self):
        with self.hold_trait_notifications():
            self.degauss.value = self.degauss_default
            self.smearing.value = self.smearing_default
            self.override.value = False


class KpointSettings(ipw.VBox):
    kpoints_distance_description = ipw.HTML(
        """<div>
        The k-points mesh density of the SCF calculation is set by the <b>protocol</b>.
        The value below represents the maximum distance between the k-points in each direction of reciprocal space.
        Tick the box to override the default, smaller is more accurate and costly. </div>"""
    )

    # The default of `kpoints_distance` must be linked to the `protocol`
    kpoints_distance_default = tl.Float(default_value=0.15)

    def __init__(self, **kwargs):
        self.override = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        self.distance = ipw.FloatText(
            value=self.kpoints_distance_default,
            step=0.05,
            description="K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override, "value"),
            (self.distance, "disabled"),
            lambda override: not override,
        )
        self.distance.observe(self.set_kpoints_distance, "value")
        self.override.observe(self.set_kpoints_distance, "value")
        self.observe(self.set_kpoints_distance, "kpoints_distance_default")

        super().__init__(
            children=[
                self.kpoints_distance_description,
                ipw.HBox([self.override, self.distance]),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    def set_kpoints_distance(self, _=None):
        self.distance.value = (
            self.distance.value
            if self.override.value
            else self.kpoints_distance_default
        )

    def _update_settings(self, **kwargs):
        """Update the kpoints_distance value by the given keyword arguments.
        This is the same as the `set_kpoints_distance` method but without the observer.
        """
        self.override.value = True
        if "kpoints_distance" in kwargs:
            self.distance.value = kwargs["kpoints_distance"]

    def reset(self):
        with self.hold_trait_notifications():
            self.distance.value = self.kpoints_distance_default
            self.override.value = False
