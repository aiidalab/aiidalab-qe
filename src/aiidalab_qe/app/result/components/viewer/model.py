from aiidalab_qe.app.result.utils import HasProcessModels, ResultsSubModel
from aiidalab_qe.common.panel import ResultsModel


class WorkflowResultsViewerModel(
    ResultsSubModel,
    HasProcessModels[ResultsModel],
):
    identifier = "results"
