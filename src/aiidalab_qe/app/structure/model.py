import traitlets as tl

from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.wizard import ConfirmableWizardStepModel, State


class StructureStepModel(
    ConfirmableWizardStepModel,
    HasInputStructure,
):
    identifier = "structure"

    structure_name = tl.Unicode("")
    manager_output = tl.Unicode("")

    installing_sssp = tl.Bool(False)
    sssp_installed = tl.Bool(allow_none=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.confirmation_exceptions += [
            "installing_sssp",
            "sssp_installed",
        ]

    def get_model_state(self) -> dict:
        return {"uuid": self.structure_uuid}

    def set_model_state(self, state: dict):
        self.structure_uuid = state.get("uuid")

    def update_state(self):
        if self.confirmed:
            self.state = State.SUCCESS
        elif self.structure_uuid:
            # We check the UUID directly (not using `has_structure`), as the structure
            # may not yet be stored.
            self.state = State.CONFIGURED
        else:
            self.state = State.READY

    def reset(self):
        self.structure_uuid = None
        self.structure_name = ""
        self.manager_output = ""

    def _check_blockers(self):
        if not self.sssp_installed:
            yield "The SSSP library is not installed"
