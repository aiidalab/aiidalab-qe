from importlib import resources

import yaml

from aiidalab_qe.common.panel import PluginOutline
from aiidalab_qe.plugins import xas as xas_folder

from .code import XasResourceSettingsModel, XasResourceSettingsPanel
from .model import XasConfigurationSettingsModel
from .result import XasResultsModel, XasResultsPanel
from .setting import XasConfigurationSettingsPanel
from .workchain import workchain_and_builder

PSEUDO_TOC = yaml.safe_load(resources.read_text(xas_folder, "pseudo_toc.yaml"))


class XasPluginOutline(PluginOutline):
    title = "X-ray absorption spectroscopy (XAS)"


xas = {
    "outline": XasPluginOutline,
    "configuration": {
        "panel": XasConfigurationSettingsPanel,
        "model": XasConfigurationSettingsModel,
    },
    "code": {
        "panel": XasResourceSettingsPanel,
        "model": XasResourceSettingsModel,
    },
    "result": {
        "panel": XasResultsPanel,
        "model": XasResultsModel,
    },
    "workchain": workchain_and_builder,
}
