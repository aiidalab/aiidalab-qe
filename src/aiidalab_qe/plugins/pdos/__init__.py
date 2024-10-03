from aiidalab_qe.common.panel import OutlinePanel
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget

from .result import Result
from .setting import Setting
from .workchain import workchain_and_builder


class PdosOutline(OutlinePanel):
    title = "Projected Density of States (PDOS)"
    help = """"""


dos_code = QEAppComputationalResourcesWidget(
    description="dos.x",
    default_calc_job_plugin="quantumespresso.dos",
)

projwfc_code = QEAppComputationalResourcesWidget(
    description="projwfc.x",
    default_calc_job_plugin="quantumespresso.projwfc",
)


pdos = {
    "outline": PdosOutline,
    "code": {"dos": dos_code, "projwfc": projwfc_code},
    "setting": Setting,
    "result": Result,
    "workchain": workchain_and_builder,
}
