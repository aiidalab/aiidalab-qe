"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.app.configuration.basic.model import WorkChainModel
from aiidalab_qe.common.panel import SettingsPanel


class WorkChainSettings(SettingsPanel[WorkChainModel]):
    title = "Basic Settings"
    identifier = "workchain"

    def __init__(self, model: WorkChainModel, **kwargs):
        super().__init__(model, **kwargs)
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

    def render(self):
        if self.rendered:
            return

        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(
            options=[("Off", "none"), ("On", "collinear")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "spin_type"),
            (self.spin_type, "value"),
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "electronic_type"),
            (self.electronic_type, "value"),
        )

        # Work chain protocol
        self.protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
        )
        ipw.link(
            (self._model, "protocol"),
            (self.protocol, "value"),
        )

        self.children = [
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                    Below you can indicate both if the material should be treated as an
                    insulator or a metal (if in doubt, choose "Metal"), and if it
                    should be studied with magnetization/spin polarization, switch
                    magnetism On or Off (On is at least twice more costly).
                </div>
            """),
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
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Protocol</h4>
                </div>
            """),
            ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
            self.protocol,
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
                    The "moderate" protocol represents a trade-off between accuracy and
                    speed. Choose the "fast" protocol for a faster calculation with
                    less precision and the "precise" protocol to aim at best accuracy
                    (at the price of longer/costlier calculations).
                </div>
            """),
        ]

        self.rendered = True

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")
