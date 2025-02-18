import traitlets as tl

from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.wizard import QeConfirmableWizardStepModel


class StructureStepModel(
    QeConfirmableWizardStepModel,
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

    def update_widget_text(self):
        if not self.has_structure:
            self.structure_name = ""
        else:
            self.manager_output = ""
            self.structure_name = str(self.input_structure.get_formula())

    def reset(self):
        self.input_structure = None
        self.structure_name = ""
        self.manager_output = ""

    def _check_blockers(self):
        if not self.sssp_installed:
            yield "The SSSP library is not installed"
