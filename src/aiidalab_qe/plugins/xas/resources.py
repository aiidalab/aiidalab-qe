"""Panel for XAS plugin."""

from aiidalab_qe.common.code.model import CodeModel, PwCodeModel
from aiidalab_qe.common.panel import ResourceSettingsModel, ResourceSettingsPanel


class XasResourceSettingsModel(ResourceSettingsModel):
    """Model for the XAS plugin."""

    identifier = "xas"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_models(
            {
                "pw": PwCodeModel(
                    name="pw.x",
                    description="pw.x",
                    default_calc_job_plugin="quantumespresso.pw",
                ),
                "xspectra": CodeModel(
                    name="xspectra.x",
                    description="xspectra.x",
                    default_calc_job_plugin="quantumespresso.xspectra",
                ),
            }
        )


class XasResourceSettingsPanel(ResourceSettingsPanel[XasResourceSettingsModel]):
    title = "XAS Structure"
