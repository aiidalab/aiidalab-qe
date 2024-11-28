from aiidalab_qe.common.panel import PluginOutline

from .model import XpsConfigurationSettingsModel
from .result import XpsResultsModel, XpsResultsPanel
from .setting import XpsConfigurationSettingsPanel
from .workchain import workchain_and_builder


class XpsPluginOutline(PluginOutline):
    title = "X-ray photoelectron spectroscopy (XPS)"


xps = {
    "outline": XpsPluginOutline,
    "configuration": {
        "panel": XpsConfigurationSettingsPanel,
        "model": XpsConfigurationSettingsModel,
    },
    "result": {
        "panel": XpsResultsPanel,
        "model": XpsResultsModel,
    },
    "workchain": workchain_and_builder,
}
