"""Panel for Bands plugin."""

from aiidalab_qe.common.code.model import CodeModel, PwCodeModel
from aiidalab_qe.common.panel import CodeSettingsModel, CodeSettingsPanel


class BandsCodeModel(CodeSettingsModel):
    """Model for the band structure plugin."""

    codes = {
        "pw": PwCodeModel(
            name="pw.x",
            description="pw.x",
            default_calc_job_plugin="quantumespresso.pw",
        ),
        "projwfc_bands": CodeModel(
            name="projwfc.x",
            description="projwfc.x",
            default_calc_job_plugin="quantumespresso.projwfc",
        ),
    }


class BandsCodeSettings(CodeSettingsPanel[BandsCodeModel]):
    title = "Band Structure"
    identifier = "bands"
