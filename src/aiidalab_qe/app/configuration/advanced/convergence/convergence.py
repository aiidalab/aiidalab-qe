import ipywidgets as ipw

from aiidalab_qe.common.widgets import HBoxWithUnits

from ..subsettings import AdvancedConfigurationSubSettingsPanel
from .model import ConvergenceConfigurationSettingsModel


class ConvergenceConfigurationSettingsPanel(
    AdvancedConfigurationSubSettingsPanel[ConvergenceConfigurationSettingsModel],
):
    def __init__(self, model: ConvergenceConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_protocol_change,
            "protocol",
        )

    def render(self):
        if self.rendered:
            return

        self.scf_conv_thr = ipw.BoundedFloatText(
            min=1e-18,
            max=1.0,
            description="Energy:",
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

        scf_conv_thr_abs = ipw.Label(
            layout=ipw.Layout(
                width="150px",
                text_align="center",
            )
        )
        ipw.dlink(
            (self._model, "scf_conv_thr"),
            (scf_conv_thr_abs, "value"),
            lambda value: f"{value * len(self._model.input_structure.sites):.5e}",
        )
        scf_conv_thr_abs.add_class("convergence-label")

        self.etot_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
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

        etot_conv_thr_abs = ipw.Label()
        ipw.dlink(
            (self._model, "etot_conv_thr"),
            (etot_conv_thr_abs, "value"),
            lambda value: f"{value * len(self._model.input_structure.sites):.5e}",
        )
        etot_conv_thr_abs.add_class("convergence-label")

        self.forc_conv_thr = ipw.BoundedFloatText(
            min=1e-15,
            max=1.0,
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

        self.mixing_mode = ipw.Dropdown(
            description="Mixing mode:",
            style={"description_width": "150px"},
        )
        ipw.dlink(
            (self._model, "mixing_mode_options"),
            (self.mixing_mode, "options"),
        )
        ipw.link(
            (self._model, "mixing_mode"),
            (self.mixing_mode, "value"),
        )

        self.mixing_beta = ipw.BoundedFloatText(
            min=0.1,
            max=1.0,
            step=0.01,
            description="Mixing beta:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "mixing_beta"),
            (self.mixing_beta, "value"),
        )

        self.help_message = ipw.HTML()
        ipw.dlink(
            (self._model, "help_message"),
            (self.help_message, "value"),
        )

        self.children = [
            self.help_message,
            ipw.HTML("<h4>Threshold for SCF cycles</h4>"),
            ipw.VBox(
                children=[
                    HBoxWithUnits(self.scf_conv_thr, "Ry/atom"),
                    HBoxWithUnits(
                        widget=scf_conv_thr_abs,
                        units="Ry",
                        layout={"margin": "-8px 0 0 150px"},
                    ),
                ]
            ),
            ipw.HTML("<h4>Thresholds for ionic convergence</h4>"),
            ipw.VBox(
                children=[
                    HBoxWithUnits(self.etot_conv_thr, "Ry/atom"),
                    HBoxWithUnits(
                        widget=etot_conv_thr_abs,
                        units="Ry",
                        layout={"margin": "-8px 0 0 150px"},
                    ),
                ]
            ),
            HBoxWithUnits(self.forc_conv_thr, "Ry/Bohr"),
            ipw.HTML("<h4>Maximum cycle steps</h4>"),
            ipw.HTML("""
                <div style="line-height: 1.4; margin-bottom: 5px;">
                    Setting a maximum number of electronic and ionic convergence steps
                    ensures that the calculation does not run indefinitely.
                </div>
            """),
            self.electron_maxstep,
            self.optimization_maxsteps,
            ipw.HTML("<h4>Mixing mode</h4>"),
            ipw.HTML("""
                <div style="line-height: 1.4; margin-bottom: 5px;">
                    The mixing mode determines how the charge density is updated during
                    the SCF cycles.
                </div>
            """),
            self.mixing_mode,
            self.mixing_beta,
        ]

        self.rendered = True

    def _on_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")
