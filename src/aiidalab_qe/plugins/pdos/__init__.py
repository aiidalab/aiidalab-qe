from aiidalab_qe.common.panel import CodePanel, OutlinePanel

from .result import Result
from .workchain import workchain_and_builder


class PdosOutline(OutlinePanel):
    title = "Projected density of states"
    help = """"""


class DosCodePanel(CodePanel):
    title = "Code"
    description = "dos.x"
    default_calc_job_plugin = "quantumespresso.dos"


class ProjwfcCodePanel(CodePanel):
    title = "Code"
    description = "projwfc.x"
    default_calc_job_plugin = "quantumespresso.projwfc"


pdos = {
    "outline": PdosOutline,
    "code": {"dos": DosCodePanel, "projwfc": ProjwfcCodePanel},
    "result": Result,
    "workchain": workchain_and_builder,
}
