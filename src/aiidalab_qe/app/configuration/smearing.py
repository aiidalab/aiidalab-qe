import ipywidgets as ipw

from .model import ConfigurationModel


class SmearingSettings(ipw.VBox):
    def __init__(self, model: ConfigurationModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[LoadingWidget("Loading smearing settings widget")],
            **kwargs,
        )

        self._model = model

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
            (self._model.advanced.smearing, "type"),
            (self.smearing, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
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
            (self._model.advanced.smearing, "degauss"),
            (self.degauss, "value"),
        )
        ipw.dlink(
            (self._model.advanced, "override"),
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

    def reset(self):
        self._model.advanced.smearing.reset()
