import ipywidgets as ipw

from aiidalab_qe.common.widgets import HBoxWithUnits

from ..subsettings import AdvancedConfigurationSubSettingsPanel
from .model import SmearingConfigurationSettingsModel


class SmearingConfigurationSettingsPanel(
    AdvancedConfigurationSubSettingsPanel[SmearingConfigurationSettingsModel],
):
    def __init__(self, model: SmearingConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_protocol_change,
            "protocol",
        )

    def render(self):
        if self.rendered:
            return

        self.smearing = ipw.Dropdown(
            description="Type:",
            style={"description_width": "150px"},
        )
        ipw.dlink(
            (self._model, "type_options"),
            (self.smearing, "options"),
        )
        ipw.link(
            (self._model, "type"),
            (self.smearing, "value"),
        )

        self.degauss = ipw.FloatText(
            step=0.005,
            description="Width:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "degauss"),
            (self.degauss, "value"),
        )

        self.children = [
            ipw.HTML("<h2>Smearing</h2>"),
            ipw.HTML("""
                <div style="line-height: 1.4; margin-bottom: 5px;">
                    Smear electronic state occupations near the Fermi level to
                    simulate finite temperature.
                    <br>
                    This helps to stabilize the SCF calculation and is important for metallic systems.
                    <br>
                    The smearing type and width are set by the chosen <b>protocol</b>.
                    <br>
                    Changes are not advised unless you've mastered
                    <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
                    target="_blank"><b>smearing effects</b></a>.
                </div>
            """),
            self.smearing,
            HBoxWithUnits(self.degauss, "Ry"),
        ]

        self.rendered = True

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")
