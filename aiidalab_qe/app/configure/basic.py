import ipywidgets as ipw

from aiidalab_qe.app.panel import Panel
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS as QEAPP_DEFAULT_PARAMETERS

DEFAULT_PARAMETERS = QEAPP_DEFAULT_PARAMETERS["basic"]


class BasicSettings(Panel):
    title = "Basic Settings"

    materials_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Below you can indicate both if the material should be treated as an insulator
        or a metal (if in doubt, choose "Metal"),
        and if it should be studied with magnetization/spin polarization,
        switch magnetism On or Off (On is at least twice more costly).
        </div>"""
    )

    protocol_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Protocol</h4></div>"""
    )
    protocol_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
        The "moderate" protocol represents a trade-off between
        accuracy and speed. Choose the "fast" protocol for a faster calculation
        with less precision and the "precise" protocol to aim at best accuracy (at the price of longer/costlier calculations).</div>"""
    )

    def __init__(self, **kwargs):
        # Work chain protocol
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value=QEAPP_DEFAULT_PARAMETERS["basic"]["protocol"],
        )
        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(
            options=[("Off", "none"), ("On", "collinear")],
            value=QEAPP_DEFAULT_PARAMETERS["basic"]["spin_type"],
            style={"description_width": "initial"},
        )
        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            value=QEAPP_DEFAULT_PARAMETERS["basic"]["electronic_type"],
            style={"description_width": "initial"},
        )

        self.children = (
            self.materials_help,
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Electronic Type:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.electronic_type,
                ]
            ),
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Magnetism:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.spin_type,
                ]
            ),
            self.protocol_title,
            ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
            self.workchain_protocol,
            self.protocol_help,
        )
        super().__init__(
            **kwargs,
        )

    def get_panel_value(self):
        """Return the value of all the widgets in the panel as a dictionary.

        :return: a dictionary of the values of all the widgets in the panel.
        """
        values = {
            "electronic_type": self.electronic_type.value,
            "spin_type": self.spin_type.value,
            "protocol": self.workchain_protocol.value,
        }
        return values

    def set_panel_value(self, parameters):
        """Load a dictionary to set the value of the widgets in the panel.

        :param panel_value: a dictionary of the values of all the widgets in the panel.
        """
        self.spin_type.value = parameters["spin_type"]
        self.electronic_type.value = parameters["electronic_type"]
        self.workchain_protocol.value = parameters["protocol"]

    def reset(self):
        self.set_panel_value(DEFAULT_PARAMETERS)
