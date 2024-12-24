"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel

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

        pseudos_model = PseudosConfigurationSettingsModel()
        self.pseudos = PseudosConfigurationSettingsPanel(model=pseudos_model)
        model.add_model("pseudos", pseudos_model)

    def render(self):
        if self.rendered:
            return

        # clean-up workchain settings
        self.clean_workdir = ipw.Checkbox(
            description="Tick to clean-up the work directory after the calculation is finished",
            indent=False,
            layout=ipw.Layout(width="fit-content"),
        )
        ipw.link(
            (self._model, "clean_workdir"),
            (self.clean_workdir, "value"),
        )
        self.reset_to_defaults_button = ipw.Button(
            description="Reset to defaults",
            button_style="warning",
            icon="undo",
            layout=ipw.Layout(width="fit-content"),
        )
        self.reset_to_defaults_button.on_click(self._on_reset_to_defaults_button_click)

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
            (self._model, "kpoints_distance"),
            (self.kpoints_distance, "value"),
        )
        ipw.dlink(
            (self._model, "input_structure"),
            (self.kpoints_distance, "disabled"),
            lambda _: not self._model.has_pbc,
        )
        self.mesh_grid = ipw.HTML()
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
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "total_charge"),
            (self.total_charge, "value"),
        )

        # Van der Waals setting widget
        self.van_der_waals = ipw.Dropdown(
            description="Van der Waals correction:",
            style={"description_width": "initial"},
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
            description="SCF conv.:",
            style={"description_width": "initial"},
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
            description="Force conv.:",
            style={"description_width": "initial"},
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
            description="Energy conv.:",
            style={"description_width": "initial"},
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
            description="Max. electron steps:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "electron_maxstep"),
            (self.electron_maxstep, "value"),
        )

        self.pseudos.render()

        self.children = [
            InAppGuide(identifier="advanced-settings"),
            ipw.HBox(
                children=[
                    self.clean_workdir,
                    self.reset_to_defaults_button,
                ],
                layout=ipw.Layout(justify_content="space-between"),
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
                layout=ipw.Layout(justify_content="space-between"),
            ),
            self.electron_maxstep,
            self.smearing,
            ipw.HTML("""
                <div>
                    The k-points mesh density of the SCF calculation is set by the
                    <b>protocol</b>. The value below represents the maximum distance
                    between the k-points in each direction of reciprocal space. Smaller
                    is more accurate and costly.
                </div>
            """),
            ipw.HBox(
                children=[
                    self.kpoints_distance,
                    self.mesh_grid,
                ]
            ),
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
