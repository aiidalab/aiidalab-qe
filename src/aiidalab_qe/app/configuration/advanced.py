# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""
import os

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from IPython.display import clear_output, display

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.panel import Panel
from aiidalab_qe.common.setup_pseudos import PseudoFamily

from .pseudos import PseudoFamilySelector, PseudoSetter


class AdvancedSettings(Panel):
    identifier = "advanced"

    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Advanced Settings</h4></div>"""
    )
    pw_adv_description = ipw.HTML(
        """Select the advanced settings for the <b>pw.x</b> code."""
    )
    kpoints_description = ipw.HTML(
        """<div>
        The k-points mesh density of the SCF calculation is set by the <b>protocol</b>.
        The value below represents the maximum distance between the k-points in each direction of reciprocal space.
        Tick the box to override the default, smaller is more accurate and costly. </div>"""
    )

    # protocol interface
    protocol = tl.Unicode(allow_none=True)
    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    spin_type = tl.Unicode()
    electronic_type = tl.Unicode()

    # output dictionary
    value = tl.Dict()

    def __init__(self, default_protocol=None, **kwargs):
        self._default_protocol = (
            default_protocol or DEFAULT_PARAMETERS["workchain"]["protocol"]
        )

        # clean-up workchain settings
        self.clean_workdir = ipw.Checkbox(
            description="",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="20px"),
        )
        self.clean_workdir_description = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Tick to clean-up the work directory after the calculation is finished.</div>"""
        )

        # Override setting widget
        self.override_prompt = ipw.HTML("<b>&nbsp;&nbsp;Override&nbsp;</b>")
        self.override = ipw.Checkbox(
            description="",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="10%"),
        )
        self.override.observe(self._override_changed, "value")

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
        self.kpoints_distance = ipw.BoundedFloatText(
            min=0.0,
            step=0.05,
            description="K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.mesh_grid = ipw.HTML()
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
        ipw.dlink(
            (self.override, "value"),
            (self.magnetization, "disabled"),
            lambda override: not override,
        )
        self.pseudo_family_selector = PseudoFamilySelector()
        self.pseudo_setter = PseudoSetter()
        ipw.dlink(
            (self.pseudo_family_selector, "value"),
            (self.pseudo_setter, "pseudo_family"),
        )
        self.kpoints_distance.observe(self._display_mesh, "value")
        self.children = [
            self.title,
            ipw.HBox(
                [self.clean_workdir, self.clean_workdir_description],
                layout=ipw.Layout(height="50px", justify_content="flex-start"),
            ),
            ipw.HBox(
                [self.pw_adv_description, self.override_widget],
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
            ipw.HBox([self.kpoints_distance, self.mesh_grid]),
            self.pseudo_family_selector,
            self.pseudo_setter,
        ]
        super().__init__(
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

        # Default settings to trigger the callback
        self.reset()

    def _override_changed(self, change):
        """Callback function to set the override value"""
        if change["new"] is False:
            # When override is set to False, reset the widget
            self.reset()

    @tl.observe("input_structure")
    def _update_input_structure(self, change):
        if self.input_structure is not None:
            self.magnetization._update_widget(change)
            self.pseudo_setter.structure = change["new"]
            self._display_mesh()
        else:
            self.magnetization.input_structure = None
            self.pseudo_setter.structure = None

    @tl.observe("electronic_type")
    def _electronic_type_changed(self, change):
        """Input electronic_type changed, update the widget values."""
        self.magnetization.electronic_type = change["new"]

    @tl.observe("protocol")
    def _protocol_changed(self, _):
        """Input protocol changed, update the widget values."""
        self._update_settings_from_protocol(self.protocol)

    def _update_settings_from_protocol(self, protocol):
        """Update the values of sub-widgets from the given protocol, this will
        trigger the callback of the sub-widget if it is exist.
        """
        self.smearing.protocol = protocol
        self.pseudo_family_selector.protocol = protocol

        parameters = PwBaseWorkChain.get_protocol_inputs(protocol)

        self.kpoints_distance.value = parameters["kpoints_distance"]

        # The pseudo_family read from the protocol (aiida-quantumespresso plugin settings)
        # we override it with the value from the pseudo_family_selector widget
        parameters["pseudo_family"] = self.pseudo_family_selector.value

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
        self.value = kwargs

    def get_panel_value(self):
        # create the the initial_magnetic_moments as None (Default)
        # XXX: start from parameters = {} and then bundle the settings by purposes (e.g. pw, bands, etc.)
        parameters = {
            "initial_magnetic_moments": None,
            "pw": {"parameters": {"SYSTEM": {}}},
            "clean_workdir": self.clean_workdir.value,
            "pseudo_family": self.pseudo_family_selector.value,
            "kpoints_distance": self.value.get("kpoints_distance"),
        }

        # Set total charge
        parameters["pw"]["parameters"]["SYSTEM"]["tot_charge"] = self.total_charge.value

        # Set the pseudos
        if self.pseudo_setter.pseudos:
            parameters["pw"]["pseudos"] = self.pseudo_setter.pseudos
            parameters["pw"]["parameters"]["SYSTEM"][
                "ecutwfc"
            ] = self.pseudo_setter.ecutwfc
            parameters["pw"]["parameters"]["SYSTEM"][
                "ecutrho"
            ] = self.pseudo_setter.ecutrho

        if self.spin_type == "collinear":
            if self.electronic_type == "metal":
                self.set_metal_parameters(parameters)
            elif self.electronic_type == "insulator":
                self.set_insulator_parameters(parameters)

        return parameters

    def set_metal_parameters(self, parameters):
        """Set the parameters for metal calculation"""
        parameters["pw"]["parameters"]["SYSTEM"][
            "smearing"
        ] = self.smearing.smearing_value
        parameters["pw"]["parameters"]["SYSTEM"][
            "degauss"
        ] = self.smearing.degauss_value
        self.set_magnetization_logic(parameters)

    def set_insulator_parameters(self, parameters):
        """Set the parameters for collinear insulator calculation"""
        parameters["pw"]["parameters"]["SYSTEM"][
            "tot_magnetization"
        ] = self.magnetization.tot_magnetization.value

    def set_magnetization_logic(self, parameters):
        """Set the parameters for magnetization calculation in metals"""
        magnetization_type = self.magnetization.magnetization_type.value
        if magnetization_type == "tot_magnetization":
            parameters["pw"]["parameters"]["SYSTEM"][
                "tot_magnetization"
            ] = self.magnetization.tot_magnetization.value
        else:
            parameters[
                "initial_magnetic_moments"
            ] = self.magnetization.get_magnetization()

    def set_panel_value(self, parameters):
        """Set the panel value from the given parameters."""

        if "pseudo_family" in parameters:
            pseudo_family_string = parameters["pseudo_family"]
            self.pseudo_family_selector.load_from_pseudo_family(
                PseudoFamily.from_string(pseudo_family_string)
            )
        if "pseudos" in parameters["pw"]:
            self.pseudo_setter.set_pseudos(parameters["pw"]["pseudos"], {})
            self.pseudo_setter.ecutwfc_setter.value = parameters["pw"]["parameters"][
                "SYSTEM"
            ]["ecutwfc"]
            self.pseudo_setter.ecutrho_setter.value = parameters["pw"]["parameters"][
                "SYSTEM"
            ]["ecutrho"]
        #
        self.kpoints_distance.value = parameters.get("kpoints_distance", 0.15)
        if parameters.get("pw") is not None:
            system = parameters["pw"]["parameters"]["SYSTEM"]
            if "degauss" in system:
                self.smearing.degauss.value = system["degauss"]
            if "smearing" in system:
                self.smearing.smearing.value = system["smearing"]
            self.total_charge.value = parameters["pw"]["parameters"]["SYSTEM"].get(
                "tot_charge", 0
            )
        if parameters.get("initial_magnetic_moments"):
            self.magnetization._set_magnetization_values(
                parameters.get("initial_magnetic_moments")
            )

    def reset(self):
        """Reset the widget and the traitlets"""

        with self.hold_trait_notifications():
            # Reset protocol dependent settings
            self._update_settings_from_protocol(self.protocol)

            # reset the pseudo family
            self.pseudo_family_selector.reset()

            # reset total charge
            self.total_charge.value = DEFAULT_PARAMETERS["advanced"]["tot_charge"]

            # reset the override checkbox
            self.override.value = False
            self.smearing.reset()
            # reset the pseudo setter
            if self.input_structure is None:
                self.pseudo_setter.structure = None
                self.pseudo_setter._reset()
            else:
                self.pseudo_setter._reset()
            # reset the magnetization
            self.magnetization.reset()
            # reset mesh grid
            if self.input_structure is None:
                self.mesh_grid.value = " "

    def _display_mesh(self, _=None):
        if self.input_structure is None:
            return
        if self.kpoints_distance.value > 0:
            # To avoid creating an aiida node every time we change the kpoints_distance,
            # we use the function itself instead of the decorated calcfunction.
            mesh = create_kpoints_from_distance.process_class._func(
                self.input_structure,
                orm.Float(self.kpoints_distance.value),
                orm.Bool(False),
            )
            self.mesh_grid.value = "Mesh " + str(mesh.get_kpoints_mesh()[0])
        else:
            self.mesh_grid.value = "Please select a number higher than 0.0"


class MagnetizationSettings(ipw.VBox):
    """Widget to set the initial magnetic moments for each kind names defined in the StructureData (StructureDtaa.get_kind_names())
    Usually these are the names of the elements in the StructureData
    (For example 'C' , 'N' , 'Fe' . However the StructureData can have defined kinds like 'Fe1' and 'Fe2')

    The widget generate a dictionary that can be used to set initial_magnetic_moments in the builder of PwBaseWorkChain

    Attributes:
        input_structure(StructureData): trait that containes the input_strucgure (confirmed structure from previous step)
    """

    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    electronic_type = tl.Unicode()
    disabled = tl.Bool()
    _DEFAULT_TOT_MAGNETIZATION = 0.0
    _DEFAULT_DESCRIPTION = "<b>Magnetization: Input structure not confirmed</b>"

    def __init__(self, **kwargs):
        self.input_structure = orm.StructureData()
        self.input_structure_labels = []
        self.tot_magnetization = ipw.BoundedFloatText(
            min=0,
            max=100,
            step=1,
            value=self._DEFAULT_TOT_MAGNETIZATION,
            disabled=True,
            description="Total magnetization:",
            style={"description_width": "initial"},
        )
        self.magnetization_type = ipw.ToggleButtons(
            options=[
                ("Tot. Magnetization", "tot_magnetization"),
                ("Starting Magnetization", "atomic_type"),
            ],
            value="tot_magnetization",
            style={"description_width": "initial"},
        )
        self.description = ipw.HTML(self._DEFAULT_DESCRIPTION)
        self.kinds = self.create_kinds_widget()
        self.kinds_widget_out = ipw.Output()
        self.magnetization_out = ipw.Output()
        self.magnetization_type.observe(self._render, "value")
        super().__init__(
            children=[
                self.description,
                self.magnetization_out,
                self.kinds_widget_out,
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    @tl.observe("disabled")
    def _disabled_changed(self, _):
        """Disable the widget"""
        if hasattr(self.kinds, "children") and self.kinds.children:
            for i in range(len(self.kinds.children)):
                self.kinds.children[i].disabled = self.disabled
        self.tot_magnetization.disabled = self.disabled
        self.magnetization_type.disabled = self.disabled

    def reset(self):
        self.disabled = True
        self.tot_magnetization.value = self._DEFAULT_TOT_MAGNETIZATION
        #
        if self.input_structure is None:
            self.description.value = self._DEFAULT_DESCRIPTION
            self.kinds = None
        else:
            self.description.value = "<b>Magnetization</b>"
            self.kinds = self.create_kinds_widget()

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

    @tl.observe("electronic_type")
    def _spin_type_changed(self, change):
        with self.magnetization_out:
            clear_output()
            if change["new"] == "metal":
                display(self.magnetization_type)
                self._render({"new": self.magnetization_type.value})
            else:
                display(self.tot_magnetization)
                with self.kinds_widget_out:
                    clear_output()

    def update_kinds_widget(self):
        self.input_structure_labels = self.input_structure.get_kind_names()
        self.kinds = self.create_kinds_widget()
        self.description.value = "<b>Magnetization</b>"

    def _render(self, value):
        if value["new"] == "tot_magnetization":
            with self.kinds_widget_out:
                clear_output()
                display(self.tot_magnetization)
        else:
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

    def _set_magnetization_values(self, magnetic_moments):
        """Set magnetization"""
        # self.override.value = True
        with self.hold_trait_notifications():
            for i in range(len(self.kinds.children)):
                if isinstance(magnetic_moments, dict):
                    self.kinds.children[i].value = magnetic_moments.get(
                        self.kinds.children[i].description, 0.0
                    )
                else:
                    self.kinds.children[i].value = magnetic_moments


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
        self._default_protocol = (
            default_protocol or DEFAULT_PARAMETERS["workchain"]["protocol"]
        )

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
            self.disabled = True
