import traitlets as tl

from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiidalab_qe.common.panel import SettingsModel


class PdosModel(SettingsModel):
    """Model for the PDOS plugin."""

    protocol = tl.Unicode(allow_none=True)

    kpoints_distance = tl.Float(0.1)
    mesh_grid = tl.Unicode("")

    def update(self):
        parameters = PdosWorkChain.get_protocol_inputs(self.protocol)
        self.kpoints_distance = parameters["nscf"]["kpoints_distance"]

    def get_model_state(self):
        return {
            "nscf_kpoints_distance": self.kpoints_distance,
        }

    def set_model_state(self, parameters: dict):
        self.nscf_kpoints_distance = parameters.get(
            "nscf_kpoints_distance",
            self.traits()["nscf_kpoints_distance"].default_value,
        )

    def reset(self):
        self.kpoints_distance = self.traits()["kpoints_distance"].default_value
        self.mesh_grid = self.traits()["mesh_grid"].default_value
