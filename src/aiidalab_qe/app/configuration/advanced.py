"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import os

import ipywidgets as ipw
import numpy as np
import traitlets as tl
from IPython.display import clear_output, display

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.panel import Panel
from aiidalab_qe.common.widgets import HubbardWidget
from aiidalab_qe.setup.pseudos import PseudoFamily

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

    dftd3_version = {
        "dft-d3": 3,
        "dft-d3bj": 4,
        "dft-d3m": 5,
        "dft-d3mbj": 6,
    }
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
        self.create_kpoints_distance_link()
        self.kpoints_distance.observe(self._callback_value_set, "value")

        # Hubbard setting widget
        self.hubbard_widget = HubbardWidget()
        ipw.dlink(
            (self.override, "value"),
            (self.hubbard_widget.activate_hubbard, "disabled"),
            lambda override: not override,
        )
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

        # Van der Waals setting widget
        self.van_der_waals = ipw.Dropdown(
            options=[
                ("None", "none"),
                ("Grimme-D3", "dft-d3"),
                ("Grimme-D3BJ", "dft-d3bj"),
                ("Grimme-D3M", "dft-d3m"),
                ("Grimme-D3MBJ", "dft-d3mbj"),
                ("Tkatchenko-Scheffler", "ts-vdw"),
            ],
            description="Van der Waals correction:",
            value="none",
            disabled=False,
            style={"description_width": "initial"},
        )

        ipw.dlink(
            (self.override, "value"),
            (self.van_der_waals, "disabled"),
            lambda override: not override,
        )

        self.magnetization = MagnetizationSettings()
        ipw.dlink(
            (self.override, "value"),
            (self.magnetization, "disabled"),
            lambda override: not override,
        )

        # Convergence Threshold settings
        self.scf_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            step=1e-10,
            description="SCF conv.:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.scf_conv_thr.observe(self._callback_value_set, "value")
        ipw.dlink(
            (self.override, "value"),
            (self.scf_conv_thr, "disabled"),
            lambda override: not override,
        )
        self.forc_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            step=0.0001,
            description="Force conv.:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.forc_conv_thr.observe(self._callback_value_set, "value")
        ipw.dlink(
            (self.override, "value"),
            (self.forc_conv_thr, "disabled"),
            lambda override: not override,
        )
        self.etot_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            step=0.00001,
            description="Energy conv.:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.etot_conv_thr.observe(self._callback_value_set, "value")
        ipw.dlink(
            (self.override, "value"),
            (self.etot_conv_thr, "disabled"),
            lambda override: not override,
        )

        # Max electron SCF steps widget
        self._create_electron_maxstep_widgets()

        # Spin-Orbit calculation
        self.spin_orbit = ipw.ToggleButtons(
            options=[
                ("Off", "wo_soc"),
                ("On", "soc"),
            ],
            description="Spin-Orbit:",
            value="wo_soc",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override, "value"),
            (self.spin_orbit, "disabled"),
            lambda override: not override,
        )

        self.pseudo_family_selector = PseudoFamilySelector()
        self.pseudo_setter = PseudoSetter()
        ipw.dlink(
            (self.pseudo_family_selector, "value"),
            (self.pseudo_setter, "pseudo_family"),
        )
        self.kpoints_distance.observe(self._display_mesh, "value")

        # Link with PseudoWidget
        ipw.dlink(
            (self.spin_orbit, "value"),
            (self.pseudo_family_selector, "spin_orbit"),
        )
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
            # van der waals setting widget
            self.van_der_waals,
            # magnetization setting widget
            self.magnetization,
            # convergence threshold setting widget
            ipw.HTML("<b>Convergence Thresholds:</b>"),
            ipw.HBox(
                [self.forc_conv_thr, self.etot_conv_thr, self.scf_conv_thr],
                layout=ipw.Layout(height="50px", justify_content="flex-start"),
            ),
            # Max electron SCF steps widget
            self.electron_maxstep,
            # smearing setting widget
            self.smearing,
            # Kpoints setting widget
            self.kpoints_description,
            ipw.HBox([self.kpoints_distance, self.mesh_grid]),
            self.hubbard_widget,
            # Spin-Orbit calculation
            self.spin_orbit,
            self.pseudo_family_selector,
            self.pseudo_setter,
        ]
        super().__init__(
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

        # Default settings to trigger the callback
        self.reset()

    def create_kpoints_distance_link(self):
        """Create the dlink for override and kpoints_distance."""
        self.kpoints_distance_link = ipw.dlink(
            (self.override, "value"),
            (self.kpoints_distance, "disabled"),
            lambda override: not override,
        )

    def remove_kpoints_distance_link(self):
        """Remove the kpoints_distance_link."""
        if hasattr(self, "kpoints_distance_link"):
            self.kpoints_distance_link.unlink()
            del self.kpoints_distance_link

    def _create_electron_maxstep_widgets(self):
        self.electron_maxstep = ipw.BoundedIntText(
            min=20,
            max=1000,
            step=1,
            value=80,
            description="Max. electron steps:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override, "value"),
            (self.electron_maxstep, "disabled"),
            lambda override: not override,
        )
        self.electron_maxstep.observe(self._callback_value_set, "value")

    def set_value_and_step(self, attribute, value):
        """
        Sets the value and adjusts the step based on the order of magnitude of the value.
        This is used for the thresolds values (etot_conv_thr, scf_conv_thr, forc_conv_thr).
        Parameters:
            attribute: The attribute whose values are to be set (e.g., self.etot_conv_thr).
            value: The numerical value to set.
        """
        attribute.value = value
        if value != 0:
            order_of_magnitude = np.floor(np.log10(abs(value)))
            attribute.step = 10 ** (order_of_magnitude - 1)
        else:
            attribute.step = 0.1  # Default step if value is zero

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
            self._update_settings_from_protocol(self.protocol)
            self._display_mesh()
            self.hubbard_widget.update_widgets(change["new"])
            if isinstance(self.input_structure, HubbardStructureData):
                self.override.value = True
            if self.input_structure.pbc == (False, False, False):
                self.kpoints_distance.value = 100.0
                self.kpoints_distance.disabled = True
                if hasattr(self, "kpoints_distance_link"):
                    self.remove_kpoints_distance_link()
            else:
                # self.kpoints_distance.disabled = False
                if not hasattr(self, "kpoints_distance_link"):
                    self.create_kpoints_distance_link()
        else:
            self.magnetization.input_structure = None
            self.pseudo_setter.structure = None
            self.hubbard_widget.update_widgets(None)
            self.kpoints_distance.disabled = False
            if not hasattr(self, "kpoints_distance_link"):
                self.create_kpoints_distance_link()

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

        if self.input_structure:
            if self.input_structure.pbc == (False, False, False):
                self.kpoints_distance.value = 100.0
                self.kpoints_distance.disabled = True
            else:
                self.kpoints_distance.value = parameters["kpoints_distance"]
        else:
            self.kpoints_distance.value = parameters["kpoints_distance"]

        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        etot_value = num_atoms * parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        self.set_value_and_step(self.etot_conv_thr, etot_value)

        # Set SCF conversion threshold
        scf_value = num_atoms * parameters["meta_parameters"]["conv_thr_per_atom"]
        self.set_value_and_step(self.scf_conv_thr, scf_value)

        # Set force conversion threshold
        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        self.set_value_and_step(self.forc_conv_thr, forc_value)

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
            "pw": {
                "parameters": {
                    "SYSTEM": {},
                    "CONTROL": {},
                    "ELECTRONS": {},
                }
            },
            "clean_workdir": self.clean_workdir.value,
            "pseudo_family": self.pseudo_family_selector.value,
            "kpoints_distance": self.value.get("kpoints_distance"),
        }

        # Set total charge
        parameters["pw"]["parameters"]["SYSTEM"]["tot_charge"] = self.total_charge.value

        if self.hubbard_widget.activate_hubbard.value:
            parameters["hubbard_parameters"] = self.hubbard_widget.hubbard_dict
            if self.hubbard_widget.eigenvalues_label.value:
                parameters["pw"]["parameters"]["SYSTEM"].update(
                    self.hubbard_widget.eigenvalues_dict
                )

        # add clean_workdir to the parameters
        parameters["clean_workdir"] = self.clean_workdir.value

        # add the pseudo_family to the parameters
        parameters["pseudo_family"] = self.pseudo_family_selector.value
        if self.pseudo_setter.pseudos:
            parameters["pw"]["pseudos"] = self.pseudo_setter.pseudos
            parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = (
                self.pseudo_setter.ecutwfc
            )
            parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"] = (
                self.pseudo_setter.ecutrho
            )

        if self.van_der_waals.value in ["none", "ts-vdw"]:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = (
                self.van_der_waals.value
            )
        else:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = "dft-d3"
            parameters["pw"]["parameters"]["SYSTEM"]["dftd3_version"] = (
                self.dftd3_version[self.van_der_waals.value]
            )

        # there are two choose, use link or parent
        if self.spin_type == "collinear":
            parameters["initial_magnetic_moments"] = (
                self.magnetization.get_magnetization()
            )
        parameters["kpoints_distance"] = self.value.get("kpoints_distance")
        if self.electronic_type == "metal":
            # smearing type setting
            parameters["pw"]["parameters"]["SYSTEM"]["smearing"] = (
                self.smearing.smearing_value
            )
            # smearing degauss setting
            parameters["pw"]["parameters"]["SYSTEM"]["degauss"] = (
                self.smearing.degauss_value
            )

        # Set tot_magnetization for collinear simulations.
        if self.spin_type == "collinear":
            # Conditions for metallic systems. Select the magnetization type and set the value if override is True
            if self.electronic_type == "metal" and self.override.value is True:
                self.set_metallic_magnetization(parameters)
            # Conditions for insulator systems. Default value is 0.0
            elif self.electronic_type == "insulator":
                self.set_insulator_magnetization(parameters)

        # convergence threshold setting
        parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"] = (
            self.forc_conv_thr.value
        )
        parameters["pw"]["parameters"]["ELECTRONS"]["conv_thr"] = (
            self.scf_conv_thr.value
        )
        parameters["pw"]["parameters"]["CONTROL"]["etot_conv_thr"] = (
            self.etot_conv_thr.value
        )

        # Max electron SCF steps
        parameters["pw"]["parameters"]["ELECTRONS"]["electron_maxstep"] = (
            self.electron_maxstep.value
        )

        # Spin-Orbit calculation
        if self.spin_orbit.value == "soc":
            parameters["pw"]["parameters"]["SYSTEM"]["lspinorb"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["noncolin"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["nspin"] = 4

        return parameters

    def set_insulator_magnetization(self, parameters):
        """Set the parameters for collinear insulator calculation. Total magnetization."""
        parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
            self.magnetization.tot_magnetization.value
        )

    def set_metallic_magnetization(self, parameters):
        """Set the parameters for magnetization calculation in metals"""
        magnetization_type = self.magnetization.magnetization_type.value
        if magnetization_type == "tot_magnetization":
            parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                self.magnetization.tot_magnetization.value
            )
        else:
            parameters["initial_magnetic_moments"] = (
                self.magnetization.get_magnetization()
            )

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
            if "lspinorb" in system:
                self.spin_orbit.value = "soc"
            else:
                self.spin_orbit.value = "wo_soc"
            # van der waals correction
            self.van_der_waals.value = self.dftd3_version.get(
                system.get("dftd3_version"),
                parameters["pw"]["parameters"]["SYSTEM"].get("vdw_corr", "none"),
            )

            # convergence threshold setting
            self.forc_conv_thr.value = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("forc_conv_thr", 0.0)
            )
            self.etot_conv_thr.value = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("etot_conv_thr", 0.0)
            )
            self.scf_conv_thr.value = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("ELECTRONS", {})
                .get("conv_thr", 0.0)
            )

            # Max electron SCF steps
            self.electron_maxstep.value = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("ELECTRONS", {})
                .get("electron_maxstep", 80)
            )

        # Logic to set the magnetization
        if parameters.get("initial_magnetic_moments"):
            self.magnetization._set_magnetization_values(
                parameters.get("initial_magnetic_moments")
            )

        if "tot_magnetization" in parameters["pw"]["parameters"]["SYSTEM"]:
            self.magnetization.magnetization_type.value = "tot_magnetization"
            self.magnetization._set_tot_magnetization(
                parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"]
            )

        if parameters.get("hubbard_parameters"):
            self.hubbard_widget.activate_hubbard.value = True
            self.hubbard_widget.set_hubbard_widget(
                parameters["hubbard_parameters"]["hubbard_u"]
            )
            starting_ns_eigenvalue = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("SYSTEM", {})
                .get("starting_ns_eigenvalue")
            )

            if starting_ns_eigenvalue is not None:
                self.hubbard_widget.eigenvalues_label.value = True
                self.hubbard_widget.set_eigenvalues_widget(starting_ns_eigenvalue)

    def reset(self):
        """Reset the widget and the traitlets"""

        with self.hold_trait_notifications():
            # Reset protocol dependent settings
            self._update_settings_from_protocol(self.protocol)

            # reset the pseudo family
            self.pseudo_family_selector.reset()

            # reset total charge
            self.total_charge.value = DEFAULT_PARAMETERS["advanced"]["tot_charge"]

            # reset the van der waals correction
            self.van_der_waals.value = DEFAULT_PARAMETERS["advanced"]["vdw_corr"]

            # reset the override checkbox
            self.override.value = False
            self.smearing.reset()
            # reset the pseudo setter
            if self.input_structure is None:
                self.pseudo_setter.structure = None
                self.pseudo_setter._reset()
            else:
                self.pseudo_setter._reset()
                if self.input_structure.pbc == (False, False, False):
                    self.kpoints_distance.value = 100.0
                    self.kpoints_distance.disabled = True

            # reset the magnetization
            self.magnetization.reset()
            # reset the hubbard widget
            self.hubbard_widget.reset()
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
    """Widget to set the type of magnetization used in the calculation:
    1) Tot_magnetization: Total majority spin charge - minority spin charge.
    2) Starting magnetization: Starting spin polarization on atomic type 'i' in a spin polarized (LSDA or noncollinear/spin-orbit) calculation.

    For Starting magnetization you can set each kind names defined in the StructureData (StructureDtaa.get_kind_names())
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
        self.tot_magnetization = ipw.BoundedIntText(
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
                ("Starting Magnetization", "starting_magnetization"),
                ("Tot. Magnetization", "tot_magnetization"),
            ],
            value="starting_magnetization",
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
    def _electronic_type_changed(self, change):
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
        self.display_kinds()

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

    def _set_tot_magnetization(self, tot_magnetization):
        """Set the total magnetization"""
        self.tot_magnetization.value = tot_magnetization


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
