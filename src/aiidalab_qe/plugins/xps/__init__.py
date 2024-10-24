from aiidalab_qe.common.panel import SettingsOutline

from .model import XpsModel
from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class XpsOutline(SettingsOutline):
    title = "X-ray photoelectron spectroscopy (XPS)"


xps = {
    "outline": XpsOutline,
    "model": XpsModel,
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
