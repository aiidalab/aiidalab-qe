from aiidalab_qe.common.panel import OutlinePanel
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget

from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class PdosOutline(OutlinePanel):
    title = "Projected density of states"
    help = """"""


pdos = {
    "outline": PdosOutline,
    "code": lambda: {
        "dos": QEAppComputationalResourcesWidget(
            description="dos.x",
            default_calc_job_plugin="quantumespresso.dos",
        ),
        "projwfc": QEAppComputationalResourcesWidget(
            description="projwfc.x",
            default_calc_job_plugin="quantumespresso.projwfc",
        ),
    },
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
