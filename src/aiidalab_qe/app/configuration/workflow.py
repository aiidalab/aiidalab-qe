"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.common.panel import Panel

from .model import config_model as model


class WorkChainSettings(Panel):
    identifier = "workchain"

    structure_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Structure</h4></div>"""
    )
    structure_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
        You have three options:<br>
        (1) Structure as is: perform a self consistent calculation using the structure provided as input.<br>
        (2) Atomic positions: perform a full relaxation of the internal atomic coordinates. <br>
        (3) Full geometry: perform a full relaxation for both the internal atomic coordinates and the cell vectors. </div>"""
    )
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
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workchain settings widget")],
            **kwargs,
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "none"),
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ],
        )
        ipw.link(
            (model.basic, "relax_type"),
            (self.relax_type, "value"),
        )

        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(
            options=[("Off", "none"), ("On", "collinear")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (model.basic, "spin_type"),
            (self.spin_type, "value"),
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (model.basic, "electronic_type"),
            (self.electronic_type, "value"),
        )

        # Work chain protocol
        self.protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
        )
        ipw.link(
            (model.basic, "protocol"),
            (self.protocol, "value"),
        )
        self.protocol.observe(self._on_protocol_change, "value")

        self.children = [
            self.structure_title,
            self.structure_help,
            self.relax_type,
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
            self.protocol,
            self.protocol_help,
        ]

        self.rendered = True

    def reset(self):
        """Reset the panel to the default value."""
        model.basic.reset()

    def get_panel_value(self):
        return {
            "protocol": model.basic.protocol,
            "relax_type": model.basic.relax_type,
            "spin_type": model.basic.spin_type,
            "electronic_type": model.basic.electronic_type,
        }

    def set_panel_value(self, parameters):
        """Update the settings based on the given dict."""
        for key in [
            "relax_type",
            "spin_type",
            "electronic_type",
        ]:
            if key in parameters:
                setattr(model.basic, key, parameters[key])
        if "protocol" in parameters:
            model.basic.protocol = parameters["protocol"]

    def _on_protocol_change(self, change):
        model.update_from_protocol(change["new"])
