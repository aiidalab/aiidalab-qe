from aiidalab_qe.common.panel import SettingsOutline

from .code import PdosCodeModel, PdosCodeSettings
from .model import PdosModel
from .result import PdosResults, PdosResultsModel
from .setting import PdosSettings
from .workchain import workchain_and_builder


class PdosOutline(SettingsOutline):
    title = "Projected Density of States (PDOS)"


pdos = {
    "outline": PdosOutline,
    "code": {
        "panel": PdosCodeSettings,
        "model": PdosCodeModel,
    },
    "setting": {
        "panel": PdosSettings,
        "model": PdosModel,
    },
    "result": {
        "panel": PdosResults,
        "model": PdosResultsModel,
    },
    "workchain": workchain_and_builder,
}
