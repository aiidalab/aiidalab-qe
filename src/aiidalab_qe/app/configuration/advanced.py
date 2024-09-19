"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.common.panel import Panel
from aiidalab_qe.setup.pseudos import PseudoFamily

from .hubbard import HubbardSettings
from .magnetization import MagnetizationSettings
from .model import config_model as model
from .pseudos import PseudoSettings
from .smearing import SmearingSettings


class AdvancedSettings(Panel):
    identifier = "advanced"

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
            (model.advanced, "clean_workdir"),
            (self.clean_workdir, "value"),
        )
        # Override setting widget
        self.override = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (model.advanced, "override"),
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
            (model.advanced, "kpoints_distance"),
            (self.kpoints_distance, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.kpoints_distance, "disabled"),
            lambda override: not override,
        )
        self.kpoints_distance.observe(
            self._on_kpoints_distance_change,
            "value",
        )
        self.mesh_grid = ipw.HTML()
        ipw.dlink(
            (model.advanced, "mesh_grid"),
            (self.mesh_grid, "value"),
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
            (model.advanced, "total_charge"),
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
            (model.advanced, "van_der_waals"),
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
            (model.advanced, "scf_conv_thr"),
            (self.scf_conv_thr, "value"),
        )
        ipw.dlink(
            (model.advanced, "scf_conv_thr_step"),
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
            (model.advanced, "forc_conv_thr"),
            (self.forc_conv_thr, "value"),
        )
        ipw.dlink(
            (model.advanced, "forc_conv_thr_step"),
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
            (model.advanced, "etot_conv_thr"),
            (self.etot_conv_thr, "value"),
        )
        ipw.dlink(
            (model.advanced, "etot_conv_thr_step"),
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
            (model.advanced, "spin_orbit"),
            (self.spin_orbit, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.spin_orbit, "disabled"),
            lambda override: not override,
        )

        self.pseudos.render()

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

        self.rendered = True

    def reset(self):
        with self.hold_trait_notifications():
            model.advanced.reset()
            self.smearing.reset()
            self.hubbard.reset()
            self.magnetization.reset()
            self.pseudos.reset()
            model.update_from_protocol()

    def get_panel_value(self):
        # create the the initial_magnetic_moments as None (Default)
        # XXX: start from parameters = {} and then bundle the settings by purposes (e.g. pw, bands, etc.)
        parameters = {
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "tot_charge": model.advanced.total_charge,
                    },
                    "CONTROL": {
                        "forc_conv_thr": model.advanced.forc_conv_thr,
                        "etot_conv_thr": model.advanced.etot_conv_thr,
                    },
                    "ELECTRONS": {
                        "conv_thr": model.advanced.scf_conv_thr,
                    },
                }
            },
            "clean_workdir": model.advanced.clean_workdir,
            "pseudo_family": model.pseudos.family,  # TODO check this
            "kpoints_distance": model.advanced.kpoints_distance,
        }

        if model.hubbard.activate:
            parameters["hubbard_parameters"] = {"hubbard_u": model.hubbard.parameters}
            if model.hubbard.eigenvalues_label:
                parameters["pw"]["parameters"]["SYSTEM"].update(
                    {"starting_ns_eigenvalue": model.hubbard.eigenvalues}
                )

        if model.pseudos.dictionary:
            parameters["pw"]["pseudos"] = model.pseudos.dictionary
            parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = model.pseudos.ecutwfc
            parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"] = model.pseudos.ecutrho

        if model.advanced.van_der_waals in ["none", "ts-vdw"]:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = (
                model.advanced.van_der_waals
            )
        else:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = "dft-d3"
            parameters["pw"]["parameters"]["SYSTEM"]["dftd3_version"] = (
                self.dftd3_version[model.advanced.van_der_waals]
            )

        # there are two choose, use link or parent
        if model.basic.spin_type == "collinear":
            parameters["initial_magnetic_moments"] = model.magnetization.moments
        if model.basic.electronic_type == "metal":
            # smearing type setting
            parameters["pw"]["parameters"]["SYSTEM"]["smearing"] = model.smearing.type
            # smearing degauss setting
            parameters["pw"]["parameters"]["SYSTEM"]["degauss"] = model.smearing.degauss

        # Set tot_magnetization for collinear simulations.
        if model.basic.spin_type == "collinear":
            # Conditions for metallic systems. Select the magnetization type and set the value if override is True
            if model.basic.electronic_type == "metal" and model.advanced.override:
                if model.magnetization.type == "tot_magnetization":
                    parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                        model.magnetization.total
                    )
                else:
                    parameters["initial_magnetic_moments"] = model.magnetization.moments
            # Conditions for insulator systems. Default value is 0.0
            elif model.basic.electronic_type == "insulator":
                parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                    model.magnetization.total
                )

        # Spin-Orbit calculation
        if model.advanced.spin_orbit == "soc":
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
        model.advanced.kpoints_distance = parameters.get("kpoints_distance", 0.15)
        if parameters.get("pw") is not None:
            system = parameters["pw"]["parameters"]["SYSTEM"]
            if "degauss" in system:
                model.smearing.degauss = system["degauss"]
            if "smearing" in system:
                model.smearing.type = system["smearing"]
            model.advanced.total_charge = parameters["pw"]["parameters"]["SYSTEM"].get(
                "tot_charge", 0
            )
            model.advanced.spin_orbit = "soc" if "lspinorb" in system else "wo_soc"
            # van der waals correction
            model.advanced.van_der_waals = self.dftd3_version.get(
                system.get("dftd3_version"),
                parameters["pw"]["parameters"]["SYSTEM"].get("vdw_corr", "none"),
            )

            # convergence threshold setting
            model.advanced.forc_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("forc_conv_thr", 0.0)
            )
            model.advanced.etot_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("etot_conv_thr", 0.0)
            )
            model.advanced.scf_conv_thr = (
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

    def _on_override_change(self, change):
        if not change["new"]:
            self.reset()

    def _on_kpoints_distance_change(self, change):
        model.advanced.update_kpoints_mesh(model.input_structure, change["new"])
