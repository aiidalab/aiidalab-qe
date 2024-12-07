"""Panel for PDOS plugin."""

from aiidalab_qe.common.code.model import CodeModel, PwCodeModel
from aiidalab_qe.common.panel import ResourceSettingsModel, ResourceSettingsPanel


class PdosResourceSettingsModel(ResourceSettingsModel):
    """Model for the pdos code setting plugin."""

    identifier = "pdos"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_models(
            {
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
        )


class PdosResourceSettingsPanel(ResourceSettingsPanel[PdosResourceSettingsModel]):
    title = "PDOS"
