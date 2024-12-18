"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.app.configuration.basic.model import BasicConfigurationSettingsModel
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel


class BasicConfigurationSettingsPanel(
    ConfigurationSettingsPanel[BasicConfigurationSettingsModel],
):
    title = "Basic Settings"
    identifier = "workchain"

    def __init__(self, model: BasicConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

    def render(self):
        if self.rendered:
            return

        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(style={"description_width": "initial"})
        ipw.dlink(
            (self._model, "spin_type_options"),
            (self.spin_type, "options"),
        )
        ipw.link(
            (self._model, "spin_type"),
            (self.spin_type, "value"),
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(style={"description_width": "initial"})
        ipw.dlink(
            (self._model, "electronic_type_options"),
            (self.electronic_type, "options"),
        )
        ipw.link(
            (self._model, "electronic_type"),
            (self.electronic_type, "value"),
        )

        # Spin-Orbit calculation
        self.spin_orbit = ipw.ToggleButtons(style={"description_width": "initial"})
        ipw.dlink(
            (self._model, "spin_orbit_options"),
            (self.spin_orbit, "options"),
        )
        ipw.link(
            (self._model, "spin_orbit"),
            (self.spin_orbit, "value"),
        )

        # Work chain protocol
        self.protocol = ipw.ToggleButtons()
        ipw.dlink(
            (self._model, "protocol_options"),
            (self.protocol, "options"),
        )
        ipw.link(
            (self._model, "protocol"),
            (self.protocol, "value"),
        )

        self.children = [
            InAppGuide(identifier="basic-settings"),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                    Below you can indicate the following:
                    <ol>
                        <li>
                            If the material should be treated as an insulator or a metal
                            (if in doubt, choose "Metal")
                        </li>
                        <li>
                            If the material should be studied with magnetization/spin
                            polarization (at least twice as costly if activated)
                        </li>
                        <li>
                            If the material should be studied with spin-orbit coupling
                        </li>
                        <li>
                            The protocol to use for the calculation, which sets default
                            values balancing the accuracy and speed of the calculation
                        </li>
                    </ol>
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
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Spin-orbit coupling:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.spin_orbit,
                ]
            ),
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Protocol:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.protocol,
                ]
            ),
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
