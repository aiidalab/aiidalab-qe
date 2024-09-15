"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiida_quantumespresso.common.types import RelaxType
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
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

    def __init__(self, callback, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workchain settings widget")],
            **kwargs,
        )

        self._set_properties(callback)

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
            (model, "relax_type"),
            (self.relax_type, "value"),
        )

        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(
            options=[("Off", "none"), ("On", "collinear")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "spin_type"),
            (self.spin_type, "value"),
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (model, "electronic_type"),
            (self.electronic_type, "value"),
        )

        # Work chain protocol
        self.protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
        )
        ipw.link(
            (model, "protocol"),
            (self.protocol, "value"),
        )

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
            *self.property_children,
            self.protocol_title,
            ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
            self.protocol,
            self.protocol_help,
        ]

        self.rendered = True

    def get_panel_value(self):
        # Work chain settings
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

        if RelaxType(model.relax_type) is not RelaxType.NONE or not (
            run_bands or run_pdos
        ):
            properties.append("relax")
        return {
            "protocol": model.protocol,
            "relax_type": model.relax_type,
            "properties": properties,
            "spin_type": model.spin_type,
            "electronic_type": model.electronic_type,
        }

    def set_panel_value(self, parameters):
        """Update the settings based on the given dict."""
        for key in [
            "relax_type",
            "spin_type",
            "electronic_type",
        ]:
            if key in parameters:
                setattr(model, key, parameters[key])
        if "protocol" in parameters:
            model.protocol = parameters["protocol"]

        properties = parameters.get("properties", [])
        for name in self.properties:
            if name in properties:
                self.properties[name].run.value = True
            else:
                self.properties[name].run.value = False

    def reset(self):
        """Reset the panel to the default value."""
        for key in ["relax_type", "spin_type", "electronic_type"]:
            setattr(model, key, model.traits()[key].default_value)
        model.protocol = model.traits()["protocol"].default_value
        for key, p in self.properties.items():
            p.run.value = key in DEFAULT_PARAMETERS["workchain"]["properties"]

    def _set_properties(self, callback):
        """Handle plugin specific settings."""

        self.properties = {}
        self.reminder_info = {}
        self.property_children = [
            self.properties_title,
            ipw.HTML("Select which properties to calculate:"),
        ]
        entries = get_entry_items("aiidalab_qe.properties", "outline")
        setting_entries = get_entry_items("aiidalab_qe.properties", "setting")
        for name, entry_point in entries.items():
            self.properties[name] = entry_point()
            self.properties[name].run.observe(callback, "value")
            self.reminder_info[name] = ipw.HTML()
            self.property_children.append(
                ipw.HBox([self.properties[name], self.reminder_info[name]])
            )

            # observer change to update the reminder text
            def update_reminder_info(change, name=name):
                if change["new"]:
                    self.reminder_info[
                        name
                    ].value = (
                        f"""Customize {name} settings in the panel above if needed."""
                    )
                else:
                    self.reminder_info[name].value = ""

            if name in setting_entries:
                self.properties[name].run.observe(update_reminder_info, "value")

        self.property_children.append(self.properties_help)
