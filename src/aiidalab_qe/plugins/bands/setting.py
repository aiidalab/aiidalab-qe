# -*- coding: utf-8 -*-
"""Panel for Bands plugin."""

import ipywidgets as ipw
import traitlets as tl
from aiidalab_qe.common.panel import Panel
from aiida import orm

class Setting(Panel):
    title = "Bands Structure"
    identifier = "bands"
    input_structure = tl.Instance(orm.StructureData, allow_none=True)

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
            description="Flat bands calculation",
            value=False,
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
        self.kpath_2d.layout.visibility = "hidden"
        self.children = [
            self.settings_title,
            self.properties_help,
            self.projwfc_bands,
            self.kpath_2d,
        ]
        super().__init__(**kwargs)


    @tl.observe("input_structure")
    def _update_structure(self, _=None):
        if self.input_structure.pbc == (True, True, False):
            self.properties_help.value = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            Please select one of the five 2D Bravais lattices corresponding to your system.
            </div>"""
            self.kpath_2d.visibility = "visible"
        elif self.input_structure.pbc == (True, False, False):
            self.properties_help.value = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            The band structure path for systems with periodicity x is from Gamma to X.
            </div>"""
            self.kpath_2d.visibility = "hidden"
        else:
            self.kpath_2d.visibility = "hidden"

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "kpath_2d": self.kpath_2d.value,
            "projwfc_bands": self.projwfc_bands.value,
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.kpath_2d.value = input_dict.get("kpath_2d", "hexagonal")
        self.projwfc_bands.value = input_dict.get("projwfc_bands", False)

    def reset(self):
        """Reset the panel to its default values."""
        self.kpath_2d.value = "hexagonal"
        self.projwfc_bands.value = False
