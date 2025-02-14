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
    def __init__(self, model: BasicConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

    def render(self):
        if self.rendered:
            return

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
        self.electronic_type.observe(
            self._on_spin_type_change,
            "value",
        )

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
        self.spin_type.observe(
            self._on_spin_type_change,
            "value",
        )

        self.starting_magmom_warning = ipw.HTML(
            value="""
                <div class="alert alert-inline alert-info" style="margin-left: 10px;">
                    <b>Note:</b> only total magnetization can be set for insulators
                </div>
            """,
            layout=ipw.Layout(display="none"),
        )

        self.magnetization_info = ipw.HTML(
            value="""
                <div class="alert alert-inline alert-info" style="margin-left: 10px;">
                    <b>Note:</b> set the desired magnetic configuration in <b>advanced
                    </b> settings
                </div>
            """,
            layout=ipw.Layout(display="none"),
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

        self.warning = ipw.HTML(
            value="""
                <div
                    class="alert alert-warning"
                    style="line-height: 140%; margin: 10px 0 0"
                >
                    <p>
                        <b>Warning:</b> detected multiples atoms with different tags.
                        You may be interested in an antiferromagnetic system. Note that
                        default starting magnetic moments do not distinguish tagged
                        atoms and are set to the same value.
                    </p>
                    <p>
                        Please go to <b>Advanced settings</b> and override the default
                        values, specifying appropriate magnetic moments for each
                        species (e.g. with different signs for an antiferromagnetic
                        configuration).
                    </p>
                </div>
            """,
            layout=ipw.Layout(display="none"),
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
                        "Electronic type:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.electronic_type,
                    self.starting_magmom_warning,
                ],
                layout=ipw.Layout(align_items="baseline"),
            ),
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Magnetism:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.spin_type,
                    self.magnetization_info,
                ],
                layout=ipw.Layout(align_items="baseline"),
            ),
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Spin-orbit coupling:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.spin_orbit,
                ],
                layout=ipw.Layout(align_items="baseline"),
            ),
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Protocol:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.protocol,
                ],
                layout=ipw.Layout(align_items="baseline"),
            ),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
                    The "moderate" protocol represents a trade-off between accuracy and
                    speed. Choose the "fast" protocol for a faster calculation with
                    less precision and the "precise" protocol to aim at best accuracy
                    (at the price of longer/costlier calculations).
                </div>
            """),
            self.warning,
        ]

        self.rendered = True

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_electronic_type_change(self, _):
        self._update_info_warning_messages()

    def _on_spin_type_change(self, _):
        self._update_info_warning_messages()

    def _update_info_warning_messages(self):
        if self._model.spin_type == "collinear":
            self.magnetization_info.layout.display = "block"
            self.starting_magmom_warning.layout.display = (
                "block" if self._model.electronic_type == "insulator" else "none"
            )
            if self._model.has_tags:
                self.warning.layout.display = "flex"
        else:
            self.starting_magmom_warning.layout.display = "none"
            self.magnetization_info.layout.display = "none"
            self.warning.layout.display = "none"
