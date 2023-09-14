# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""
import ipywidgets as ipw
from aiida_quantumespresso.common.types import RelaxType

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items


class WorkChainSettings(ipw.VBox):
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

    properties_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Properties</h4></div>"""
    )
    properties_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
        The band structure workflow will
        automatically detect the default path in reciprocal space using the
        <a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">
        SeeK-path tool</a>.</div>"""
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
        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "none"),
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ],
            value="positions_cell",
        )

        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(
            options=[("Off", "none"), ("On", "collinear")],
            value=DEFAULT_PARAMETERS["spin_type"],
            style={"description_width": "initial"},
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            value=DEFAULT_PARAMETERS["electronic_type"],
            style={"description_width": "initial"},
        )

        # Work chain protocol
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        self.properties = {}
        self.property_children = [
            self.properties_title,
            ipw.HTML("Select which properties to calculate:"),
        ]
        entries = get_entry_items("aiidalab_qe.properties", "outline")
        for name, entry_point in entries.items():
            self.properties[name] = entry_point()
            self.property_children.append(self.properties[name])
        self.property_children.append(self.properties_help)
        super().__init__(
            children=[
                self.structure_title,
                self.structure_help,
                self.relax_type,
                self.materials_help,
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Electronic Type:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="120px"
                            ),
                        ),
                        self.electronic_type,
                    ]
                ),
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Magnetism:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="120px"
                            ),
                        ),
                        self.spin_type,
                    ]
                ),
                *self.property_children,
                self.protocol_title,
                ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
                self.workchain_protocol,
                self.protocol_help,
            ],
            **kwargs,
        )

    def get_setting_parameters(self):
        # Work chain settings
        relax_type = self.relax_type.value
        electronic_type = self.electronic_type.value
        spin_type = self.spin_type.value

        protocol = self.workchain_protocol.value

        properties = []

        # add plugin specific settings
        run_bands = False
        run_pdos = False
        for name in self.properties:
            if self.properties[name].run.value:
                properties.append(name)
            if name == "bands":
                run_bands = True
            elif name == "pdos":
                run_bands = True

        if RelaxType(relax_type) is not RelaxType.NONE or not (run_bands or run_pdos):
            properties.append("relax")
        return {
            "protocol": protocol,
            "relax_type": relax_type,
            "properties": properties,
            "spin_type": spin_type,
            "electronic_type": electronic_type,
        }

    def set_setting_parameters(self, parameters):
        """Update the settings based on the given dict."""
        for key in [
            "relax_type",
            "spin_type",
            "electronic_type",
            "workchain_protocol",
        ]:
            if key in parameters:
                getattr(self, key).value = parameters[key]
        properties = parameters.get("properties", [])
        for name in self.properties:
            if name in properties:
                self.properties[name].run.value = True
            else:
                self.properties[name].run.value = False
