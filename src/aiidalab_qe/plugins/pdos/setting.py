"""Panel for Pdos plugin."""

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiidalab_qe.common.panel import Panel

RYDBERG_TO_EV = 13.605703976


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
        self.description = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
                By default, the tetrahedron method is used for PDOS calculation. If required you can apply Gaussian broadening with a custom degauss value.
                <br>
                For molecules and systems with localized orbitals, it is recommended to use a custom degauss value.
            </div>"""
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
        self.use_pdos_degauss = ipw.Checkbox(
            value=False,
            description="Use custom PDOS degauss",
            style={"description_width": "initial"},
        )
        self.pdos_degauss = ipw.BoundedFloatText(
            min=0.00001,
            step=0.001,
            value=0.005,
            description="PDOS degauss (Ry):",
            disabled=True,
            style={"description_width": "initial"},
        )
        self.pdos_degauss_eV = ipw.HTML(
            value=f"({self.pdos_degauss.value * RYDBERG_TO_EV:.4f} eV)",
        )

        self.use_pdos_degauss.observe(self._disable_pdos_degauss, "value")
        self.pdos_degauss.observe(self._update_pdos_degauss_ev, "value")
        self.mesh_grid = ipw.HTML()
        self.nscf_kpoints_distance.observe(self._display_mesh, "value")
        self.nscf_kpoints_distance.observe(self._procotol_changed, "change")
        self.children = [
            self.settings_title,
            self.description,
            ipw.HBox([self.nscf_kpoints_distance, self.mesh_grid]),
            self.use_pdos_degauss,
            ipw.HBox([self.pdos_degauss, self.pdos_degauss_eV]),
        ]
        super().__init__(**kwargs)

    def _disable_pdos_degauss(self, change):
        self.pdos_degauss.disabled = not change["new"]

    def _update_pdos_degauss_ev(self, change):
        new_value = change["new"] * RYDBERG_TO_EV
        if self.pdos_degauss_eV.value != f"({new_value} eV)":
            self.pdos_degauss_eV.value = f"({new_value:.4f} eV)"

    @tl.observe("protocol")
    def _procotol_changed(self, change):
        self.nscf_kpoints_distance.value = PdosWorkChain.get_protocol_inputs(
            change["new"]
        )["nscf"]["kpoints_distance"]
        self._display_mesh()

    @tl.observe("input_structure")
    def _update_structure(self, _=None):
        self._display_mesh()
        # For molecules this is compulsory
        if self.input_structure and self.input_structure.pbc == (False, False, False):
            self.nscf_kpoints_distance.value = 100
            self.nscf_kpoints_distance.disabled = True
            self.use_pdos_degauss.value = True
            self.use_pdos_degauss.disabled = True

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
            "pdos_degauss": self.pdos_degauss.value,
            "use_pdos_degauss": self.use_pdos_degauss.value,
        }

    def set_panel_value(self, input_dict):
        """Load a dictionary with the input parameters for the plugin."""
        self.nscf_kpoints_distance.value = input_dict.get("nscf_kpoints_distance", 0.1)
        self.pdos_degauss.value = input_dict.get("pdos_degauss", 0.005)
        self.use_pdos_degauss.value = input_dict.get("use_pdos_degauss", False)

    def reset(self):
        """Reset the panel to its default values."""
        self.nscf_kpoints_distance.value = 0.1
        self.pdos_degauss.value = 0.005
        self.use_pdos_degauss.value = False
