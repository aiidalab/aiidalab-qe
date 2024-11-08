from aiidalab_qe.common.panel import SettingsOutline

from .model import XpsModel
from .result import Result
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
    "result": Result,
    "workchain": workchain_and_builder,
}
