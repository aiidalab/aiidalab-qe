import ipywidgets as ipw

from aiidalab_qe.common.mixins import HasModels, HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.panel import ResultsModel


class WorkChainViewerModel(
    Model,
    HasModels[ResultsModel],
    HasProcess,
):
    def _link_model(self, model: ResultsModel):
        ipw.dlink(
            (self, "monitor_counter"),
            (model, "monitor_counter"),
        )
