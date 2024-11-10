import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.mixins import HasModels, HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.panel import ResultsModel


# TODO remove if structure viewer receives the MVC treatment
class StructureViewerModel(ResultsModel):
    """Stand-in model for the structure viewer."""

    include = True


class WorkChainViewerModel(
    Model,
    HasModels[ResultsModel],
    HasProcess,
):
    include = tl.Bool(False)

    def _link_model(self, model: ResultsModel):
        ipw.dlink(
            (self, "process_node"),
            (model, "process_node"),
        )
