"""Panel for Bands plugin."""

import ipywidgets as ipw

from aiidalab_qe.common.panel import SettingPanel


class Setting(SettingPanel):
    title = "Bands Structure"
    identifier = "bands"

    def render(self):
        if self.rendered:
            return

        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Settings</h4></div>"""
        )
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        self.kpath_2d_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            If your system has periodicity xy. Please select one of the five 2D Bravais lattices corresponding to your system.
            </div>"""
        )
        self.kpath_2d = ipw.Dropdown(
            description="Lattice:",
            options=[
                ("Hexagonal", "hexagonal"),
                ("Square", "square"),
                ("Rectangular", "rectangular"),
                ("Centered Rectangular", "centered_rectangular"),
                ("Oblique", "oblique"),
            ],
            value="hexagonal",
        )

        self.children = [
            self.settings_title,
            self.kpath_2d_help,
            self.kpath_2d,
        ]

        self.rendered = True

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "kpath_2d": self.kpath_2d.value,
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.kpath_2d.value = input_dict.get("kpath_2d", "hexagonal")

    def reset(self):
        """Reset the panel to its default values."""
        self.kpath_2d.value = "hexagonal"
