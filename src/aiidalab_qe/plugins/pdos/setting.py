"""Panel for Pdos plugin."""

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiidalab_qe.common.panel import Panel


class Setting(Panel):
    title = "PDOS"
    identifier = "pdos"
    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    protocol = tl.Unicode(allow_none=True)

    def __init__(self, **kwargs):
        self.settings_title = ipw.HTML(
            """<div style="padding-top: 0px; padding-bottom: 0px">
            <h4>Settings</h4></div>"""
        )
        # nscf kpoints setting widget
        self.nscf_kpoints_distance = ipw.BoundedFloatText(
            min=0.001,
            step=0.01,
            value=0.1,
            description="NSCF K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.mesh_grid = ipw.HTML()
        self.nscf_kpoints_distance.observe(self._display_mesh, "value")
        self.nscf_kpoints_distance.observe(self._procotol_changed, "change")
        self.children = [
            self.settings_title,
            ipw.HBox([self.nscf_kpoints_distance, self.mesh_grid]),
        ]
        super().__init__(**kwargs)

    @tl.observe("protocol")
    def _procotol_changed(self, change):
        self.nscf_kpoints_distance.value = PdosWorkChain.get_protocol_inputs(
            change["new"]
        )["nscf"]["kpoints_distance"]
        self._display_mesh()

    @tl.observe("input_structure")
    def _update_structure(self, _=None):
        self._display_mesh()

    def _display_mesh(self, _=None):
        if self.input_structure is None:
            return
        mesh = create_kpoints_from_distance.process_class._func(
            self.input_structure,
            orm.Float(self.nscf_kpoints_distance.value),
            orm.Bool(False),
        )
        self.mesh_grid.value = "Mesh " + str(mesh.get_kpoints_mesh()[0])

    def get_panel_value(self):
        """Return a dictionary with the input parameters for the plugin."""
        return {
            "nscf_kpoints_distance": self.nscf_kpoints_distance.value,
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.nscf_kpoints_distance.value = input_dict.get("nscf_kpoints_distance", 0.1)

    def reset(self):
        """Reset the panel to its default values."""
        self.nscf_kpoints_distance.value = 0.1
