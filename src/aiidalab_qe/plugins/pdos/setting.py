"""Panel for Pdos plugin."""

import ipywidgets as ipw

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiidalab_qe.app.configuration.model import ConfigurationModel
from aiidalab_qe.common.panel import SettingsPanel


class Setting(SettingsPanel):
    title = "PDOS"
    identifier = "pdos"

    def __init__(self, config_model: ConfigurationModel, **kwargs):
        super().__init__(config_model, **kwargs)
        ipw.dlink(
            (self._config_model.workchain, "protocol"),
            (self._model, "protocol"),
        )
        self._config_model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._config_model.workchain.observe(
            self._on_protocol_change,
            "protocol",
        )
        self._model.observe(
            self._on_kpoints_distance_change,
            "kpoints_distance",
        )

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

        self.rendered = True

    def update(self):
        if self.updated:
            return
        self._update_mesh()
        self.updated = True

    def reset(self):
        self._model.reset()
        self.updated = False

    def _on_input_structure_change(self, _=None):
        self._update_mesh()

    def _on_protocol_change(self, _):
        self._model.update()
        self._update_mesh()

    def _on_kpoints_distance_change(self, _=None):
        self._update_mesh()

    def _update_mesh(self, _=None):
        if not self._model.include or self._config_model.input_structure is None:
            self._model.mesh_grid = ""
        else:
            mesh = create_kpoints_from_distance.process_class._func(
                self._config_model.input_structure,
                orm.Float(self._model.kpoints_distance),
                orm.Bool(False),
            )
            self._model.mesh_grid = "Mesh " + str(mesh.get_kpoints_mesh()[0])
