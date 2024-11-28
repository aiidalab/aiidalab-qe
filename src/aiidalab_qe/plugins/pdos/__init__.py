from aiidalab_qe.common.panel import PluginOutline

from .code import PdosResourceSettingsModel, PdosResourceSettingsPanel
from .model import PdosConfigurationSettingsModel
from .result import PdosResultsModel, PdosResultsPanel
from .setting import PdosConfigurationSettingPanel
from .workchain import workchain_and_builder


class PdosPluginOutline(PluginOutline):
    title = "Projected Density of States (PDOS)"


pdos = {
    "outline": PdosPluginOutline,
    "configuration": {
        "panel": PdosConfigurationSettingPanel,
        "model": PdosConfigurationSettingsModel,
    },
    "code": {
        "panel": PdosResourceSettingsPanel,
        "model": PdosResourceSettingsModel,
    },
    "result": {
        "panel": PdosResultsPanel,
        "model": PdosResultsModel,
    },
    "workchain": workchain_and_builder,
}
