# from aiidalab_qe.bands.result import Result
from aiidalab_qe.common.panel import SettingsOutline

from .code import BandsCodeModel, BandsCodeSettings
from .model import BandsModel
from .result import BandsResults, BandsResultsModel
from .setting import BandsSettings
from .workchain import workchain_and_builder


class BandsOutline(SettingsOutline):
    title = "Electronic band structure"


bands = {
    "outline": BandsOutline,
    "code": {
        "panel": BandsCodeSettings,
        "model": BandsCodeModel,
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
