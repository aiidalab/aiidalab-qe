import ipywidgets as ipw
import traitlets as tl

from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.common.widgets import LoadingWidget

from .model import config_model as model


class SmearingSettings(ipw.VBox):
    smearing_description = ipw.HTML(
        """<p>
        The smearing type and width is set by the chosen <b>protocol</b>.
        Tick the box to override the default, not advised unless you've mastered <b>smearing effects</b> (click <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
        target="_blank">here</a> for a discussion).
    </p>"""
    )

    protocol = tl.Unicode(allow_none=True)

    def __init__(self, **kwargs):
        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[LoadingWidget("Loading smearing settings widget")],
            **kwargs,
        )

        self.rendered = False

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
            (model, "smearing"),
            (self.smearing, "value"),
        )
        ipw.dlink(
            (model, "override"),
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
            (model, "degauss"),
            (self.degauss, "value"),
        )
        ipw.dlink(
            (model, "override"),
            (self.degauss, "disabled"),
            lambda override: not override,
        )

        self.children = [
            self.smearing_description,
            ipw.HBox([self.smearing, self.degauss]),
        ]

        ipw.dlink(
            (model, "protocol"),
            (self, "protocol"),
        )

        self.rendered = True

    @tl.observe("protocol")
    def _update_settings_from_protocol(self, _=None):
        """Update the widget values from the given protocol, and trigger the callback."""
        parameters = PwBaseWorkChain.get_protocol_inputs(model.protocol)["pw"][
            "parameters"
        ]["SYSTEM"]
        model.degauss = parameters["degauss"]
        model.smearing = parameters["smearing"]

    def reset(self):
        """Reset the widget and the traitlets"""
        model.protocol = model.traits()["protocol"].default_value
