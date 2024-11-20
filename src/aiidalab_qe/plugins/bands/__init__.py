# from aiidalab_qe.bands.result import Result
from aiidalab_qe.app.submission.code.model import CodeModel
from aiidalab_qe.common.panel import SettingsOutline

from .model import BandsModel
from .result import BandsResults, BandsResultsModel
from .setting import BandsSettings
from .workchain import workchain_and_builder


class BandsOutline(SettingsOutline):
    title = "Electronic band structure"


bands = {
    "outline": BandsOutline,
    "code": {
        "projwfc_bands": CodeModel(
            description="projwfc.x",
            default_calc_job_plugin="quantumespresso.projwfc",
        ),
    },
    "setting": {
        "panel": BandsSettings,
        "model": BandsModel,
    },
    "result": {
        "panel": BandsResults,
        "model": BandsResultsModel,
    },
    "workchain": workchain_and_builder,
}
