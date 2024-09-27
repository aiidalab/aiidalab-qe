"""Panel for Pdos plugin."""

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiidalab_qe.common.panel import SettingPanel


class Setting(SettingPanel):
    title = "PDOS"
    identifier = "pdos"

    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    protocol = tl.Unicode(allow_none=True)

    def render(self):
        if self.rendered:
            return

        self.nscf_kpoints_distance = ipw.BoundedFloatText(
            min=0.001,
            step=0.01,
            description="NSCF K-points distance (1/Å):",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "kpoints_distance"),
            (self.nscf_kpoints_distance, "value"),
        )
        self.nscf_kpoints_distance.observe(
            self._on_kpoints_distance_change,
            "value",
        )

        self.mesh_grid = ipw.HTML()
        ipw.dlink(
            (self._model, "mesh_grid"),
            (self.mesh_grid, "value"),
        )

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Settings</h4>
                </div>
            """),
            ipw.HBox(
                children=[
                    self.nscf_kpoints_distance,
                    self.mesh_grid,
                ]
            ),
        ]

        with self.hold_trait_notifications():
            ipw.dlink(
                (self._config_model, "input_structure"),
                (self, "input_structure"),
            )
            ipw.dlink(
                (self._config_model.workchain, "protocol"),
                (self._model, "protocol"),
            )
            ipw.dlink(
                (self._model, "protocol"),
                (self, "protocol"),
            )

        self.rendered = True

    def reset(self):
        self._model.reset()

    @tl.observe("input_structure")
    def _on_input_structure_change(self, _=None):
        self._update_mesh()

    @tl.observe("protocol")
    def _on_protocol_change(self, _):
        self._model.update()
        self._update_mesh()

    def _on_kpoints_distance_change(self, _=None):
        self._update_mesh()

    def _update_mesh(self, _=None):
        if self.input_structure is None:
            self._model.mesh_grid = ""
        else:
            mesh = create_kpoints_from_distance.process_class._func(
                self.input_structure,
                orm.Float(self._model.kpoints_distance),
                orm.Bool(False),
            )
            self._model.mesh_grid = "Mesh " + str(mesh.get_kpoints_mesh()[0])
