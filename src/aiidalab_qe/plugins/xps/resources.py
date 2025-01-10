"""Resource panel for XPS plugin."""

from aiidalab_qe.common.code.model import PwCodeModel
from aiidalab_qe.common.panel import (
    PluginResourceSettingsModel,
    PluginResourceSettingsPanel,
)


class XpsResourceSettingsModel(PluginResourceSettingsModel):
    """Model for the XPS plugin."""

    title = "XPS"
    identifier = "xps"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_models(
            {
                "pw": PwCodeModel(
                    name="pw.x",
                    description="pw.x",
                    default_calc_job_plugin="quantumespresso.pw",
                ),
            }
        )


class XpsResourceSettingsPanel(
    PluginResourceSettingsPanel[XpsResourceSettingsModel],
):
    """Panel for configuring the XPS plugin."""
