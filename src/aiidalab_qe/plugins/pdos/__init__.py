from aiidalab_qe.common.panel import PluginOutline

from .model import PdosConfigurationSettingsModel
from .resources import PdosResourceSettingsModel, PdosResourceSettingsPanel
from .result import PdosResultsModel, PdosResultsPanel
from .setting import PdosConfigurationSettingPanel
from .workchain import workchain_and_builder


class PdosPluginOutline(PluginOutline):
    title = "Electronic projected density of states (PDOS)"


pdos = {
    "outline": PdosPluginOutline,
    "configuration": {
        "panel": PdosConfigurationSettingPanel,
        "model": PdosConfigurationSettingsModel,
    },
    "resources": {
        "panel": PdosResourceSettingsPanel,
        "model": PdosResourceSettingsModel,
    },
    "result": {
        "panel": PdosResultsPanel,
        "model": PdosResultsModel,
    },
    "workchain": workchain_and_builder,
}
