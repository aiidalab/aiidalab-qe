"""Panel for Bands plugin."""

import ipywidgets as ipw

from aiidalab_qe.common.panel import Panel


class Setting(Panel):
    title = "Bands Structure"
    identifier = "bands"

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Settings</h4></div>"""
        )
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        self.properties_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
        The band structure workflow will
        automatically detect the default path in reciprocal space using the
        <a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">
        SeeK-path tool</a>.</div>"""
        )
        self.projwfc_bands = ipw.Checkbox(
            description="Fat bands calculation",
            value=False,
        )
        self.children = [
            self.settings_title,
            self.properties_help,
            self.projwfc_bands,
        ]
        super().__init__(**kwargs)

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "projwfc_bands": self.projwfc_bands.value,
        }

    def set_panel_value(self, input_dict):
        self.projwfc_bands.value = input_dict.get("projwfc_bands", False)

    def reset(self):
        """Reset the panel to its default values."""
        self.projwfc_bands.value = False
