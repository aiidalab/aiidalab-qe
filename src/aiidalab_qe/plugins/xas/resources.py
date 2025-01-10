"""Panel for XAS plugin."""

from aiidalab_qe.common.code.model import CodeModel, PwCodeModel
from aiidalab_qe.common.panel import (
    PluginResourceSettingsModel,
    PluginResourceSettingsPanel,
)


class XasResourceSettingsModel(PluginResourceSettingsModel):
    """Model for the XAS plugin."""

    title = "XAS"
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


class XasResourceSettingsPanel(
    PluginResourceSettingsPanel[XasResourceSettingsModel],
):
    """Panel for configuring the XAS plugin."""
