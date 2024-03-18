from aiidalab_qe.common.panel import OutlinePanel

from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class XpsOutline(OutlinePanel):
    title = "X-ray photoelectron spectroscopy (XPS)"
    help = """"""


xps = {
    "outline": XpsOutline,
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
