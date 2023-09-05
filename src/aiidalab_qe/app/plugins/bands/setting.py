# -*- coding: utf-8 -*-
"""Panel for Bands plugin."""
import ipywidgets as ipw
from aiida.orm import Int, Str

from aiidalab_qe.app.common.panel import Panel


class Setting(Panel):
    title = "Bands Structure Settings"
    identifier = "bands"

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Settings</h4></div>"""
        )
        self.settings_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Please set the value of path and number of points.
            </div>"""
        )
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        self.path = ipw.Text(
            value="G",
            description="Bands path:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.npoint = ipw.IntText(
            value=50,
            description="Number of point:",
            disabled=False,
            style={"description_width": "initial"},
        )

        self.widgets = [
            self.settings_title,
            self.settings_help,
            self.path,
            self.npoint,
        ]
        super().__init__(**kwargs)

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "path": Str(self.path.value),
            "npoint": Int(self.npoint.value),
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.path.value = input_dict.get("path", 1)
        self.npoint.value = input_dict.get("npoint", 2)
