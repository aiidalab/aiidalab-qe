from aiidalab_qe.app.panel import OutlinePanel

from .result import Result
from .workchain import workchain_and_builder


class PDOSOutline(OutlinePanel):
    title = "Projected density of states (PDOS)"


property = {
    "outline": PDOSOutline,
    "result": Result,
    "workchain": workchain_and_builder,
}
