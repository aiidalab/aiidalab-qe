from aiidalab_qe.common.panel import SettingsOutline

from .model import XpsModel
from .result import XpsResults, XpsResultsModel
from .setting import XpsSettings
from .workchain import workchain_and_builder


class XpsOutline(SettingsOutline):
    title = "X-ray photoelectron spectroscopy (XPS)"


xps = {
    "outline": XpsOutline,
    "setting": {
        "panel": XpsSettings,
        "model": XpsModel,
    },
    "result": {
        "panel": XpsResults,
        "model": XpsResultsModel,
    },
    "workchain": workchain_and_builder,
}
