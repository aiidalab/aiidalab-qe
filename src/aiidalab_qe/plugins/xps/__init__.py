from aiidalab_qe.common.panel import PanelOutline

from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class XpsOutline(PanelOutline):
    title = "X-ray photoelectron spectroscopy (XPS)"
    help = """"""


xps = {
    "outline": XpsOutline,
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
