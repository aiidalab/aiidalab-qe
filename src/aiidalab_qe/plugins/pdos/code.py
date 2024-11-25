"""Panel for PDOS plugin."""

from aiidalab_qe.common.code.model import CodeModel, PwCodeModel
from aiidalab_qe.common.panel import CodeSettingsModel, CodeSettingsPanel


class PdosCodeModel(CodeSettingsModel):
    """Model for the pdos code setting plugin."""

    codes = {
        "pw": PwCodeModel(
            name="pw.x",
            description="pw.x",
            default_calc_job_plugin="quantumespresso.pw",
        ),
        "dos": CodeModel(
            name="dos.x",
            description="dos.x",
            default_calc_job_plugin="quantumespresso.dos",
        ),
        "projwfc": CodeModel(
            name="projwfc.x",
            description="projwfc.x",
            default_calc_job_plugin="quantumespresso.projwfc",
        ),
    }


class PdosCodeSettings(CodeSettingsPanel[PdosCodeModel]):
    title = "PDOS"
    identifier = "pdos"
