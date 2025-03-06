"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel
from aiidalab_qe.common.widgets import HBoxWithUnits

from .hubbard import (
    HubbardConfigurationSettingsModel,
    HubbardConfigurationSettingsPanel,
)
from .magnetization import (
    MagnetizationConfigurationSettingsModel,
    MagnetizationConfigurationSettingsPanel,
)
from .model import AdvancedConfigurationSettingsModel
from .pseudos import (
    PseudosConfigurationSettingsModel,
    PseudosConfigurationSettingsPanel,
)
from .smearing import (
    SmearingConfigurationSettingsModel,
    SmearingConfigurationSettingsPanel,
)


class AdvancedConfigurationSettingsPanel(
    ConfigurationSettingsPanel[AdvancedConfigurationSettingsModel],
):
    def __init__(self, model: AdvancedConfigurationSettingsModel, **kwargs):
        super().__init__(
            model=model,
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            **kwargs,
        )

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_protocol_change,
            "protocol",
        )
        self._model.observe(
            self._on_kpoints_distance_change,
            "kpoints_distance",
        )

        # NOTE connect pseudos first, as some settings depend on it
        pseudos_model = PseudosConfigurationSettingsModel()
        self.pseudos = PseudosConfigurationSettingsPanel(model=pseudos_model)
        model.add_model("pseudos", pseudos_model)

        smearing_model = SmearingConfigurationSettingsModel()
        self.smearing = SmearingConfigurationSettingsPanel(model=smearing_model)
        model.add_model("smearing", smearing_model)

        magnetization_model = MagnetizationConfigurationSettingsModel()
        self.magnetization = MagnetizationConfigurationSettingsPanel(
            model=magnetization_model,
        )
        model.add_model("magnetization", magnetization_model)

        hubbard_model = HubbardConfigurationSettingsModel()
        self.hubbard = HubbardConfigurationSettingsPanel(model=hubbard_model)
        model.add_model("hubbard", hubbard_model)

    def render(self):
        if self.rendered:
            return

        self.reset_to_defaults_button = ipw.Button(
            description="Reset to defaults",
            button_style="primary",
            icon="undo",
            layout=ipw.Layout(width="fit-content"),
        )
        self.reset_to_defaults_button.on_click(self._on_reset_to_defaults_button_click)

        self.clean_workdir = ipw.Checkbox(
            description="Delete the work directory after the calculation",
            indent=False,
            layout=ipw.Layout(width="fit-content", margin="5px 2px"),
        )
        ipw.link(
            (self._model, "clean_workdir"),
            (self.clean_workdir, "value"),
        )

        # Smearing setting widget
        self.smearing.render()

        # Kpoints setting widget
        self.kpoints_distance = ipw.BoundedFloatText(
            min=0.0,
            step=0.05,
            description="K-points distance:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "kpoints_distance"),
            (self.kpoints_distance, "value"),
        )
        ipw.dlink(
            (self._model, "input_structure"),
            (self.kpoints_distance, "disabled"),
            lambda _: not self._model.has_pbc,
        )
        self.mesh_grid = ipw.HTML(layout=ipw.Layout(margin="0 0 0 10px"))
        ipw.dlink(
            (self._model, "mesh_grid"),
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
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "total_charge"),
            (self.total_charge, "value"),
        )

        # Van der Waals setting widget
        self.van_der_waals = ipw.Dropdown(
            description="Van der Waals correction:",
            style={"description_width": "150px"},
        )
        ipw.dlink(
            (self._model, "van_der_waals_options"),
            (self.van_der_waals, "options"),
        )
        ipw.link(
            (self._model, "van_der_waals"),
            (self.van_der_waals, "value"),
        )

        # Magnetization settings
        self.magnetization.render()

        # Convergence Threshold settings
        self.scf_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            description="SCF:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "scf_conv_thr"),
            (self.scf_conv_thr, "value"),
        )
        ipw.dlink(
            (self._model, "scf_conv_thr_step"),
            (self.scf_conv_thr, "step"),
        )
        self.forc_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            format="0.0e",
            description="Force:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "forc_conv_thr"),
            (self.forc_conv_thr, "value"),
        )
        ipw.dlink(
            (self._model, "forc_conv_thr_step"),
            (self.forc_conv_thr, "step"),
        )
        self.etot_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
            format="0.0e",
            description="Energy:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "etot_conv_thr"),
            (self.etot_conv_thr, "value"),
        )
        ipw.dlink(
            (self._model, "etot_conv_thr_step"),
            (self.etot_conv_thr, "step"),
        )
        self.electron_maxstep = ipw.BoundedIntText(
            min=20,
            max=1000,
            step=1,
            description="Electronic:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "electron_maxstep"),
            (self.electron_maxstep, "value"),
        )

        self.optimization_maxsteps = ipw.BoundedIntText(
            min=50,
            max=1000,
            step=1,
            description="Ionic:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "optimization_maxsteps"),
            (self.optimization_maxsteps, "value"),
        )
        self.pseudos.render()

        self.children = [
            InAppGuide(identifier="advanced-settings"),
            self.reset_to_defaults_button,
            self.clean_workdir,
            self.total_charge,
            self.van_der_waals,
            ipw.HTML("<h2>Convergence</h2>"),
            ipw.HTML("""
                <div style="line-height: 1.4; margin-bottom: 5px;">
                    Control the convergence criteria of the self-consistent field (SCF)
                    geometry optimization cycles.
                </div>
            """),
            ipw.HTML("<h4>Thresholds</h4>"),
            ipw.HTML("""
                <div style="line-height: 1.4; margin-bottom: 5px;">
                    Setting thresholds for energy, force, and self-consistency ensures calculation accuracy and stability.
                    <br>
                    Lower values increase the accuracy but also the computational cost.
                    <br>
                    The default values are set by the <b>protocol</b> are usually a
                    good starting point.
                </div>
            """),
            HBoxWithUnits(self.forc_conv_thr, "Ry/Bohr"),
            HBoxWithUnits(self.etot_conv_thr, "Ry/atom"),
            HBoxWithUnits(self.scf_conv_thr, "Ry/atom"),
            ipw.HTML("<h4>Maximum cycle steps</h4>"),
            ipw.HTML("""
                <div style="line-height: 1.4; margin-bottom: 5px;">
                    Setting a maximum number of electronic and ionic optimization steps
                    ensures that the calculation does not run indefinitely.
                </div>
            """),
            self.electron_maxstep,
            self.optimization_maxsteps,
            self.smearing,
            ipw.HTML("<h2>K-points</h2>"),
            ipw.HTML("""
                <div style="line-height: 1.4; margin-bottom: 5px;">
                    The k-points mesh density of the SCF calculation is set by the
                    <b>protocol</b>.
                    <br>
                    The value below represents the maximum distance between k-points
                    in each direction of reciprocal space.
                    <br>
                    Smaller is more accurate and costly.
                </div>
            """),
            ipw.HBox(
                children=[
                    HBoxWithUnits(self.kpoints_distance, "Å<sup>-1</sup>"),
                    self.mesh_grid,
                ],
                layout=ipw.Layout(align_items="center"),
            ),
            self.magnetization,
            self.hubbard,
            self.pseudos,
        ]

        self.rendered = True

        self.refresh()

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")

    def _on_kpoints_distance_change(self, _):
        self.refresh(specific="mesh")

    def _on_reset_to_defaults_button_click(self, _):
        self._reset()

    def _reset(self):
        self._model.reset()
        for _, model in self._model.get_models():
            model.reset()
