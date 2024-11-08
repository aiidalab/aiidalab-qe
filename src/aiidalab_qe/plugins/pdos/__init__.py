from aiidalab_qe.app.submission.code import CodeModel
from aiidalab_qe.common.panel import SettingsOutline

from .model import PdosModel
from .result import Result
from .setting import PdosSettings
from .workchain import workchain_and_builder


class PdosOutline(SettingsOutline):
    title = "Projected Density of States (PDOS)"


pdos = {
    "outline": PdosOutline,
    "code": {
        "dos": CodeModel(
            description="dos.x",
            default_calc_job_plugin="quantumespresso.dos",
        ),
        "projwfc": CodeModel(
            description="projwfc.x",
            default_calc_job_plugin="quantumespresso.projwfc",
        ),
    },
    "setting": {
        "panel": PdosSettings,
        "model": PdosModel,
    },
    "result": Result,
    "workchain": workchain_and_builder,
}
