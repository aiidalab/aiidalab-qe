import traitlets as tl

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import PanelModel


class PdosConfigurationSettingsModel(PanelModel, HasInputStructure):
    title = "PDOS"
    identifier = "pdos"

    dependencies = [
        "structure_uuid",
        "workchain.protocol",
    ]

    protocol = tl.Unicode(allow_none=True)

    nscf_kpoints_distance = tl.Float(0.1)
    mesh_grid = tl.Unicode("")
    use_pdos_degauss = tl.Bool(False)
    pdos_degauss = tl.Float(0.005)
    energy_grid_step = tl.Float(0.01)

    def update(self, specific=""):
        with self.hold_trait_notifications():
            if not specific or specific != "mesh":
                parameters = PdosWorkChain.get_protocol_inputs(self.protocol)
                self._update_nscf_kpoints_distance(parameters)

            self._update_kpoints_mesh()

    def get_model_state(self) -> dict:
        return {
            "nscf_kpoints_distance": self.nscf_kpoints_distance,
            "use_pdos_degauss": self.use_pdos_degauss,
            "pdos_degauss": self.pdos_degauss,
            "energy_grid_step": self.energy_grid_step,
        }

    def set_model_state(self, state: dict):
        self.nscf_kpoints_distance = state.get(
            "nscf_kpoints_distance",
            self._get_default("nscf_kpoints_distance"),
        )
        self.use_pdos_degauss = state.get("use_pdos_degauss", False)
        self.pdos_degauss = state.get("pdos_degauss", 0.005)
        self.energy_grid_step = state.get("energy_grid_step", 0.01)

    def reset(self):
        with self.hold_trait_notifications():
            self.nscf_kpoints_distance = self._get_default("nscf_kpoints_distance")
            self.use_pdos_degauss = self._get_default("use_pdos_degauss")
            self.pdos_degauss = self._get_default("pdos_degauss")
            self.energy_grid_step = self._get_default("energy_grid_step")

    def _update_kpoints_mesh(self, _=None):
        if not self.has_structure:
            mesh_grid = ""
        elif self.nscf_kpoints_distance > 0:
            mesh = create_kpoints_from_distance.process_class._func(
                self.input_structure,
                orm.Float(self.nscf_kpoints_distance),
                orm.Bool(False),
            )
            mesh_grid = f"Mesh {mesh.get_kpoints_mesh()[0]!s}"
        else:
            mesh_grid = "Please select a number higher than 0.0"
        self._defaults["mesh_grid"] = mesh_grid
        self.mesh_grid = mesh_grid

    def _update_nscf_kpoints_distance(self, parameters):
        if self.has_pbc:
            nscf_kpoints_distance = parameters["nscf"]["kpoints_distance"]
        else:
            nscf_kpoints_distance = 100.0
            self.use_pdos_degauss = True
        self._defaults["nscf_kpoints_distance"] = nscf_kpoints_distance
        self.nscf_kpoints_distance = self._defaults["nscf_kpoints_distance"]
