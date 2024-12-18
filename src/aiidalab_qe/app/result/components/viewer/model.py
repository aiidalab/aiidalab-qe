import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponentModel
from aiidalab_qe.common.mixins import HasModels
from aiidalab_qe.common.panel import ResultsModel


class WorkChainResultsViewerModel(
    ResultsComponentModel,
    HasModels[ResultsModel],
):
    identifier = "workflow results"

    def _link_model(self, model: ResultsModel):
        ipw.dlink(
            (self, "process_uuid"),
            (model, "process_uuid"),
        )
        ipw.dlink(
            (self, "monitor_counter"),
            (model, "monitor_counter"),
        )
