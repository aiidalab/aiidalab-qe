from aiidalab_qe.app.result.utils import HasProcessModels, ResultsSubModel
from aiidalab_qe.common.mixins import HasProcess


class WorkChainStatusModel(
    ResultsSubModel,
    HasProcessModels[HasProcess],
):
    identifier = "status"
