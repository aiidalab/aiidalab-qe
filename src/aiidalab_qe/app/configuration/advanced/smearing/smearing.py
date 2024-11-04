import ipywidgets as ipw

from ..subsettings import AdvancedSubSettings
from .model import SmearingModel


class SmearingSettings(AdvancedSubSettings[SmearingModel]):
    identifier = "smearing"

    def __init__(self, model: SmearingModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_protocol_change,
            "protocol",
        )

    def render(self):
        if self.rendered:
            return

        self.smearing = ipw.Dropdown(
            options=["cold", "gaussian", "fermi-dirac", "methfessel-paxton"],
            description="Smearing type:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "type"),
            (self.smearing, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.smearing, "disabled"),
            lambda override: not override,
        )

        self.degauss = ipw.FloatText(
            step=0.005,
            description="Smearing width (Ry):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "degauss"),
            (self.degauss, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.degauss, "disabled"),
            lambda override: not override,
        )

        self.children = [
            ipw.HTML("""
                <p>
                    The smearing type and width is set by the chosen <b>protocol</b>.
                    Tick the box to override the default, not advised unless you've
                    mastered <b>smearing effects</b> (click
                    <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
                    target="_blank">here</a> for a discussion).
                </p>
            """),
            ipw.HBox(
                children=[
                    self.smearing,
                    self.degauss,
                ]
            ),
        ]

        self.rendered = True

        self.refresh()

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")
