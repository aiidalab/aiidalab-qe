import ipywidgets as ipw
import traitlets as tl

from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS


class SmearingSettings(ipw.VBox):
    # accept protocol as input and set the values
    protocol = tl.Unicode(allow_none=True)

    # The output of the widget is a dictionary with the values of smearing and degauss
    degauss_value = tl.Float()
    smearing_value = tl.Unicode()

    smearing_description = ipw.HTML(
        """<p>
        The smearing type and width is set by the chosen <b>protocol</b>.
        Tick the box to override the default, not advised unless you've mastered <b>smearing effects</b> (click <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
        target="_blank">here</a> for a discussion).
    </p>"""
    )
    disabled = tl.Bool()

    def __init__(self, default_protocol=None, **kwargs):
        self._default_protocol = (
            default_protocol or DEFAULT_PARAMETERS["workchain"]["protocol"]
        )

        self.smearing = ipw.Dropdown(
            options=["cold", "gaussian", "fermi-dirac", "methfessel-paxton"],
            description="Smearing type:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.degauss = ipw.FloatText(
            step=0.005,
            description="Smearing width (Ry):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self, "disabled"),
            (self.degauss, "disabled"),
        )
        ipw.dlink(
            (self, "disabled"),
            (self.smearing, "disabled"),
        )
        self.degauss.observe(self._callback_value_set, "value")
        self.smearing.observe(self._callback_value_set, "value")

        super().__init__(
            children=[
                self.smearing_description,
                ipw.HBox([self.smearing, self.degauss]),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

        # Default settings to trigger the callback
        self.protocol = self._default_protocol

    @tl.default("disabled")
    def _default_disabled(self):
        return False

    @tl.observe("protocol")
    def _protocol_changed(self, _):
        """Input protocol changed, update the widget values."""
        self._update_settings_from_protocol(self.protocol)

    def _update_settings_from_protocol(self, protocol):
        """Update the widget values from the given protocol, and trigger the callback."""
        parameters = PwBaseWorkChain.get_protocol_inputs(protocol)["pw"]["parameters"][
            "SYSTEM"
        ]

        with self.hold_trait_notifications():
            # This changes will trigger callbacks
            self.degauss.value = parameters["degauss"]
            self.smearing.value = parameters["smearing"]

    def _callback_value_set(self, _=None):
        """callback function to set the smearing and degauss values"""
        settings = {
            "degauss": self.degauss.value,
            "smearing": self.smearing.value,
        }
        self.update_settings(**settings)

    def update_settings(self, **kwargs):
        """Set the output dict from the given keyword arguments.
        This function will only update the traitlets but not the widget value.
        """
        self.degauss_value = kwargs.get("degauss")
        self.smearing_value = kwargs.get("smearing")

    def reset(self):
        """Reset the widget and the traitlets"""
        self.protocol = self._default_protocol

        with self.hold_trait_notifications():
            self._update_settings_from_protocol(self.protocol)
            self.disabled = True
