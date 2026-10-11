import traitlets as tl

from aiidalab_qe.app.result.utils import HasProcessModels, ResultsSubModel
from aiidalab_qe.common.panel import ResultsModel


class WorkflowResultsViewerModel(
    ResultsSubModel,
    HasProcessModels[ResultsModel],
):
    identifier = "results"
    plugin_status_checked = tl.Bool(False)
    incompatible_plugin_ids = tl.Set(tl.Unicode())
    incompatible_plugin_titles = tl.List(tl.Unicode())
