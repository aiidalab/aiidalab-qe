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
from aiidalab_qe.common.panel import Panel
from aiidalab_qe.setup.pseudos import PseudoFamily

from .hubbard import HubbardSettings
from .magnetization import MagnetizationSettings
from .model import config_model as model
from .pseudos import PseudoSettings
from .smearing import SmearingSettings


class AdvancedSettings(Panel):
    identifier = "advanced"

    protocol = tl.Unicode(allow_none=True)
    input_structure = tl.Instance(orm.StructureData, allow_none=True)

    def __init__(self, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[LoadingWidget("Loading advanced settings widget")],
            **kwargs,
        )

        self.smearing = SmearingSettings()
        self.hubbard = HubbardSettings()
        self.magnetization = MagnetizationSettings()
        self.pseudos = PseudoSettings()

        self.dftd3_version = {
            "dft-d3": 3,
            "dft-d3bj": 4,
            "dft-d3m": 5,
            "dft-d3mbj": 6,
        }

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        # clean-up workchain settings
        self.clean_workdir = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="20px"),
        )
        ipw.link(
            (model, "clean_workdir"),
            (self.clean_workdir, "value"),
        )
        # Override setting widget
        self.override = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (model, "override"),
            (self.override, "value"),
        )
        ipw.dlink(
            (model, "input_structure"),
            (self.override, "disabled"),
            lambda structure: structure is None,
        )
        self.override.observe(self._on_override_change, "value")

        # Smearing setting widget
        self.smearing.render()

        # Kpoints setting widget
        self.kpoints_distance = ipw.BoundedFloatText(
            min=0.0,
            step=0.05,
            description="K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "kpoints_distance"),
            (self.kpoints_distance, "value"),
        )
        self.mesh_grid = ipw.HTML()
        ipw.dlink(
            (self.override, "value"),
            (self.kpoints_distance, "disabled"),
            lambda override: not override,
        )

        # Hubbard setting widget
        self.hubbard.render()

        # Total change setting widget
        self.total_charge = ipw.BoundedFloatText(
            min=-3,
            max=3,
            step=0.01,
            disabled=False,
            description="Total charge:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "total_charge"),
            (self.total_charge, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.total_charge, "disabled"),
            lambda override: not override,
        )

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
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "van_der_waals"),
            (self.van_der_waals, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.van_der_waals, "disabled"),
            lambda override: not override,
        )

        # Magnetization settings
        self.magnetization.render()

        # Convergence Threshold settings
        self.scf_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            description="SCF conv.:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "scf_conv_thr"),
            (self.scf_conv_thr, "value"),
        )
        ipw.dlink(
            (model, "scf_conv_thr_step"),
            (self.scf_conv_thr, "step"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.scf_conv_thr, "disabled"),
            lambda override: not override,
        )
        self.forc_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            description="Force conv.:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "forc_conv_thr"),
            (self.forc_conv_thr, "value"),
        )
        ipw.dlink(
            (model, "forc_conv_thr_step"),
            (self.forc_conv_thr, "step"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.forc_conv_thr, "disabled"),
            lambda override: not override,
        )
        self.etot_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            description="Energy conv.:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "etot_conv_thr"),
            (self.etot_conv_thr, "value"),
        )
        ipw.dlink(
            (model, "etot_conv_thr_step"),
            (self.etot_conv_thr, "step"),
        )
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
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "spin_orbit"),
            (self.spin_orbit, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.spin_orbit, "disabled"),
            lambda override: not override,
        )

        self.pseudos.render()

        self.kpoints_distance.observe(self._display_mesh, "value")

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 10px">
                    <h4>Advanced Settings</h4>
                </div>
            """),
            ipw.HBox(
                children=[
                    self.clean_workdir,
                    ipw.HTML("""
                        <div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
                            Tick to clean-up the work directory after the calculation is finished.
                        </div>
                    """),
                ],
                layout=ipw.Layout(height="50px", justify_content="flex-start"),
            ),
            ipw.HBox(
                children=[
                    ipw.HTML("""
                        Select the advanced settings for the <b>pw.x</b> code.
                    """),
                    ipw.HBox(
                        children=[
                            ipw.HTML(
                                value="<b>Override</b>",
                                layout=ipw.Layout(margin="0 5px 0 0"),
                            ),
                            self.override,
                        ],
                        layout=ipw.Layout(max_width="20%"),
                    ),
                ],
                layout=ipw.Layout(height="50px", justify_content="space-between"),
            ),
            self.total_charge,
            self.van_der_waals,
            self.magnetization,
            ipw.HTML("<b>Convergence Thresholds:</b>"),
            ipw.HBox(
                children=[
                    self.forc_conv_thr,
                    self.etot_conv_thr,
                    self.scf_conv_thr,
                ],
                layout=ipw.Layout(height="50px", justify_content="flex-start"),
            ),
            self.smearing,
            ipw.HTML("""
                <div>
                    The k-points mesh density of the SCF calculation is set by the
                    <b>protocol</b>. The value below represents the maximum distance
                    between the k-points in each direction of reciprocal space. Tick
                    the box to override the default, smaller is more accurate and
                    costly.
                </div>
            """),
            ipw.HBox(
                children=[
                    self.kpoints_distance,
                    self.mesh_grid,
                ]
            ),
            self.hubbard,
            self.spin_orbit,
            self.pseudos,
        ]

        ipw.dlink(
            (model, "input_structure"),
            (self, "input_structure"),
        )
        ipw.dlink(
            (model, "protocol"),
            (self, "protocol"),
        )

        self.rendered = True

    def reset(self):
        """Reset the widget and the traitlets"""

        with self.hold_trait_notifications():
            model.protocol = model.traits()["protocol"].default_value
            model.total_charge = model.traits()["total_charge"].default_value
            model.van_der_waals = model.traits()["van_der_waals"].default_value
            model.override = False
            self.smearing.reset()
            self.hubbard.reset()
            self.magnetization.reset()
            self.pseudos.reset()
            if model.input_structure is None:
                model.mesh_grid = " "

    def get_panel_value(self):
        # create the the initial_magnetic_moments as None (Default)
        # XXX: start from parameters = {} and then bundle the settings by purposes (e.g. pw, bands, etc.)
        parameters = {
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "tot_charge": model.total_charge,
                    },
                    "CONTROL": {
                        "forc_conv_thr": model.forc_conv_thr,
                        "etot_conv_thr": model.etot_conv_thr,
                    },
                    "ELECTRONS": {
                        "conv_thr": model.scf_conv_thr,
                    },
                }
            },
            "clean_workdir": model.clean_workdir,
            "pseudo_family": model.pseudos.family,  # TODO check this
            "kpoints_distance": model.kpoints_distance,
        }

        if model.hubbard.activate:
            parameters["hubbard_parameters"] = {"hubbard_u": model.hubbard.parameters}
            if model.hubbard.eigenvalues_label:
                parameters["pw"]["parameters"]["SYSTEM"].update(
                    {"starting_ns_eigenvalue": model.hubbard.eigenvalues}
                )

        if model.pseudos.dict:
            parameters["pw"]["pseudos"] = model.pseudos.dict
            parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = model.pseudos.ecutwfc
            parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"] = model.pseudos.ecutrho

        if model.van_der_waals in ["none", "ts-vdw"]:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = model.van_der_waals
        else:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = "dft-d3"
            parameters["pw"]["parameters"]["SYSTEM"]["dftd3_version"] = (
                self.dftd3_version[model.van_der_waals]
            )

        # there are two choose, use link or parent
        if model.spin_type == "collinear":
            parameters["initial_magnetic_moments"] = model.magnetization.moments
        if model.electronic_type == "metal":
            # smearing type setting
            parameters["pw"]["parameters"]["SYSTEM"]["smearing"] = model.smearing.type
            # smearing degauss setting
            parameters["pw"]["parameters"]["SYSTEM"]["degauss"] = model.smearing.degauss

        # Set tot_magnetization for collinear simulations.
        if model.spin_type == "collinear":
            # Conditions for metallic systems. Select the magnetization type and set the value if override is True
            if model.electronic_type == "metal" and model.override:
                if model.magnetization.type == "tot_magnetization":
                    parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                        model.magnetization.total
                    )
                else:
                    parameters["initial_magnetic_moments"] = model.magnetization.moments
            # Conditions for insulator systems. Default value is 0.0
            elif model.electronic_type == "insulator":
                parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                    model.magnetization.total
                )

        # Spin-Orbit calculation
        if model.spin_orbit == "soc":
            parameters["pw"]["parameters"]["SYSTEM"]["lspinorb"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["noncolin"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["nspin"] = 4

        return parameters

    def set_panel_value(self, parameters):
        """Set the panel value from the given parameters."""

        if "pseudo_family" in parameters:
            pseudo_family = PseudoFamily.from_string(parameters["pseudo_family"])
            library = pseudo_family.library
            accuracy = pseudo_family.accuracy
            model.pseudos.library = f"{library} {accuracy}"
            model.pseudos.functional = pseudo_family.functional
        if "pseudos" in parameters["pw"]:
            model.pseudos.dict = parameters["pw"]["pseudos"]
            model.pseudos.ecutwfc = parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
            model.pseudos.ecutrho = parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"]
        #
        model.kpoints_distance = parameters.get("kpoints_distance", 0.15)
        if parameters.get("pw") is not None:
            system = parameters["pw"]["parameters"]["SYSTEM"]
            if "degauss" in system:
                model.smearing.degauss = system["degauss"]
            if "smearing" in system:
                model.smearing.type = system["smearing"]
            model.total_charge = parameters["pw"]["parameters"]["SYSTEM"].get(
                "tot_charge", 0
            )
            model.spin_orbit = "soc" if "lspinorb" in system else "wo_soc"
            # van der waals correction
            model.van_der_waals = self.dftd3_version.get(
                system.get("dftd3_version"),
                parameters["pw"]["parameters"]["SYSTEM"].get("vdw_corr", "none"),
            )

            # convergence threshold setting
            model.forc_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("forc_conv_thr", 0.0)
            )
            model.etot_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("etot_conv_thr", 0.0)
            )
            model.scf_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("ELECTRONS", {})
                .get("conv_thr", 0.0)
            )

        # Logic to set the magnetization
        if magnetic_moments := parameters.get("initial_magnetic_moments"):
            if isinstance(magnetic_moments, list):
                magnetic_moments = {
                    kind: magnetic_moments[i]
                    for i, kind in enumerate(model.input_structure.get_kind_names())
                }
            model.magnetization.moments = magnetic_moments

        if "tot_magnetization" in parameters["pw"]["parameters"]["SYSTEM"]:
            model.magnetization.type = "tot_magnetization"

        if parameters.get("hubbard_parameters"):
            model.hubbard.activate = True
            model.hubbard.parameters = parameters["hubbard_parameters"]["hubbard_u"]
            starting_ns_eigenvalue = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("SYSTEM", {})
                .get("starting_ns_eigenvalue")
            )
            if starting_ns_eigenvalue is not None:
                model.hubbard.eigenvalues_label = True
                model.hubbard.eigenvalues = starting_ns_eigenvalue

    @tl.observe("protocol")
    def _on_protocol_change(self, _=None):
        self._update_model()

    @tl.observe("input_structure")
    def _on_input_structure_change(self, _):
        if model.input_structure:
            self._update_model()
            self._display_mesh()
            if isinstance(model.input_structure, HubbardStructureData):
                model.override = True

    def _on_override_change(self, change):
        if not change["new"]:
            self.reset()

    def _update_model(self):
        """Update the model from the given protocol."""

        def set_value_and_step(attribute, value):
            """Sets the value and step size.

            Parameters:
                attribute (str):
                    The attribute whose values are to be set (e.g., self.etot_conv_thr).
                value (float):
                    The numerical value to set.
            """
            setattr(model, attribute, value)
            if value != 0:
                order_of_magnitude = np.floor(np.log10(abs(value)))
                setattr(model, f"{attribute}_step", 10 ** (order_of_magnitude - 1))
            else:
                setattr(model, f"{attribute}_step", 0.1)

        parameters = PwBaseWorkChain.get_protocol_inputs(model.protocol)

        model.kpoints_distance = parameters["kpoints_distance"]

        num_atoms = len(model.input_structure.sites) if model.input_structure else 1

        etot_value = num_atoms * parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        set_value_and_step("etot_conv_thr", etot_value)

        scf_value = num_atoms * parameters["meta_parameters"]["conv_thr_per_atom"]
        set_value_and_step("scf_conv_thr", scf_value)

        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        set_value_and_step("forc_conv_thr", forc_value)

    def _display_mesh(self, _=None):
        if model.input_structure is None:
            return
        if model.kpoints_distance > 0:
            # To avoid creating an aiida node every time we change the kpoints_distance,
            # we use the function itself instead of the decorated calcfunction.
            mesh = create_kpoints_from_distance.process_class._func(
                model.input_structure,
                orm.Float(model.kpoints_distance),
                orm.Bool(False),
            )
            model.mesh_grid = "Mesh " + str(mesh.get_kpoints_mesh()[0])
        else:
            model.mesh_grid = "Please select a number higher than 0.0"
