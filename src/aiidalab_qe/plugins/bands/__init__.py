# from aiidalab_qe.bands.result import Result
from pathlib import Path

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
    "guides": Path(__file__).parent / "guides",
}
