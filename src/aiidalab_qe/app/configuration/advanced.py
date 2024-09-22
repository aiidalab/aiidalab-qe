"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.common.panel import Panel

from .hubbard import HubbardSettings
from .magnetization import MagnetizationSettings
from .model import ConfigurationModel
from .pseudos import PseudoSettings
from .smearing import SmearingSettings


class AdvancedSettings(Panel):
    identifier = "advanced"

    def __init__(self, model: ConfigurationModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[LoadingWidget("Loading advanced settings widget")],
            **kwargs,
        )

        self._model = model

        self.smearing = SmearingSettings(model=model)
        self.hubbard = HubbardSettings(model=model)
        self.magnetization = MagnetizationSettings(model=model)
        self.pseudos = PseudoSettings(model=model)

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
            (self._model.advanced, "clean_workdir"),
            (self.clean_workdir, "value"),
        )
        # Override setting widget
        self.override = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (self._model.advanced, "override"),
            (self.override, "value"),
        )
        ipw.dlink(
            (self._model, "input_structure"),
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
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model.advanced, "kpoints_distance"),
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
            (self._model.advanced, "mesh_grid"),
            (self.mesh_grid, "value"),
        )

        # Hubbard setting widget
        self.hubbard.render()

        # Total change setting widget
        self.total_charge = ipw.BoundedFloatText(
            min=-3,
            max=3,
            step=0.01,
            description="Total charge:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model.advanced, "total_charge"),
            (self.total_charge, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
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
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model.advanced, "van_der_waals"),
            (self.van_der_waals, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
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
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model.advanced, "scf_conv_thr"),
            (self.scf_conv_thr, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "scf_conv_thr_step"),
            (self.scf_conv_thr, "step"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
            (self.scf_conv_thr, "disabled"),
            lambda override: not override,
        )
        self.forc_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            description="Force conv.:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model.advanced, "forc_conv_thr"),
            (self.forc_conv_thr, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "forc_conv_thr_step"),
            (self.forc_conv_thr, "step"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
            (self.forc_conv_thr, "disabled"),
            lambda override: not override,
        )
        self.etot_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            description="Energy conv.:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model.advanced, "etot_conv_thr"),
            (self.etot_conv_thr, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "etot_conv_thr_step"),
            (self.etot_conv_thr, "step"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
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
            (self._model.advanced, "spin_orbit"),
            (self.spin_orbit, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
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
            self._model.advanced.reset()
            self.smearing.reset()
            self.hubbard.reset()
            self.magnetization.reset()
            self.pseudos.reset()

    def _on_kpoints_distance_change(self, _=None):
        self._model.advanced.update_kpoints_mesh()

    def _on_override_change(self, change):
        if not change["new"]:
            self.reset()
