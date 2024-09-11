"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw
import numpy as np
import traitlets as tl

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.panel import Panel
from aiidalab_qe.setup.pseudos import PseudoFamily

from .hubbard import HubbardWidget
from .magnetization import MagnetizationSettings
from .pseudos import PseudoFamilySelector, PseudoSetter
from .smearing import SmearingSettings


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
        from aiidalab_qe.common.widgets import LoadingWidget

        self._default_protocol = (
            default_protocol or DEFAULT_PARAMETERS["workchain"]["protocol"]
        )

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[LoadingWidget("Loading advanced settings widget")],
            **kwargs,
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        from .model import config_model

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

        ipw.dlink(
            (config_model, "protocol"),
            (self, "protocol"),
        )
        ipw.dlink(
            (config_model, "spin_type"),
            (self, "spin_type"),
        )
        ipw.dlink(
            (config_model, "electronic_type"),
            (self, "electronic_type"),
        )
        ipw.dlink(
            (config_model, "input_structure"),
            (self, "input_structure"),
        )

        # Default settings to trigger the callback
        self.reset()

        self.rendered = True

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
        else:
            self.magnetization.input_structure = None
            self.pseudo_setter.structure = None
            self.hubbard_widget.update_widgets(None)

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
