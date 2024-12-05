# from aiidalab_qe.bands.result import Result
from aiidalab_qe.common.panel import PluginOutline

from .model import BandsConfigurationSettingsModel
from .resources import BandsResourceSettingsModel, BandsResourceSettingsPanel
from .result import BandsResultsModel, BandsResultsPanel
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
    "result": {
        "panel": BandsResultsPanel,
        "model": BandsResultsModel,
    },
    "workchain": workchain_and_builder,
}
