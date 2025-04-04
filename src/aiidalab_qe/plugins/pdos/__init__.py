from aiidalab_qe.common.panel import PluginOutline

from .model import PdosConfigurationSettingsModel
from .resources import PdosResourceSettingsModel, PdosResourceSettingsPanel
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
    "workchain": workchain_and_builder,
    "metadata": {
        "process_labels": {
            "PdosWorkChain": "Projected density of states workflow",
            "DosCalculation": "Density of states calculation",
            "ProjwfcCalculation": "Orbital projection calculation",
        }
    },
}
