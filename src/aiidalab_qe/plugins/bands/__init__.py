from aiidalab_qe.common.panel import PluginOutline

from .model import BandsConfigurationSettingsModel
from .resources import BandsResourceSettingsModel, BandsResourceSettingsPanel
from .setting import BandsConfigurationSettingsPanel
from .workchain import workchain_and_builder


class BandsPluginOutline(PluginOutline):
    title = "Electronic band structure"


bands = {
    "outline": BandsPluginOutline,
    "configuration": {
        "panel": BandsConfigurationSettingsPanel,
        "model": BandsConfigurationSettingsModel,
    },
    "resources": {
        "panel": BandsResourceSettingsPanel,
        "model": BandsResourceSettingsModel,
    },
    "workchain": workchain_and_builder,
    "metadata": {
        "process_labels": {
            "BandsWorkChain": "Electronic band structure workflow",
            "PwBandsWorkChain": "Electronic band structure workflow",
            "ProjwfcBandsWorkChain": "Electronic band structure workflow",
            "ProjwfcCalculation": "Orbital projection calculation",
        }
    },
}
