from aiidalab_qe.common.panel import OutlinePanel

from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class PdosOutline(OutlinePanel):
    title = "Projected density of states"
    help = """"""


pdos = {
    "outline": PdosOutline,
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
