from importlib import resources

import yaml

from aiidalab_qe.common.panel import SettingsOutline
from aiidalab_qe.plugins import xas as xas_folder

from .model import XasModel
from .result import XasResults, XasResultsModel
from .setting import XasSettings
from .workchain import workchain_and_builder
from .code import XasCodeModel, XasCodeSettings

PSEUDO_TOC = yaml.safe_load(resources.read_text(xas_folder, "pseudo_toc.yaml"))


class XasOutline(SettingsOutline):
    title = "X-ray absorption spectroscopy (XAS)"


xas = {
    "outline": XasOutline,
    "code": {
        "panel": XasCodeSettings,
        "model": XasCodeModel,
    },
    "setting": {
        "panel": XasSettings,
        "model": XasModel,
    },
    "result": {
        "panel": XasResults,
        "model": XasResultsModel,
    },
    "workchain": workchain_and_builder,
}
