from aiidalab_qe.common.panel import OutlinePanel

from .result import Result
from .workchain import workchain_and_builder


class PdosOutline(OutlinePanel):
    title = "Projected density of states"
    help = """"""


pdos = {
    "outline": PdosOutline,
    "result": Result,
    "workchain": workchain_and_builder,
}
