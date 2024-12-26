"""PDOS results view widgets"""

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import PdosResultsModel


class PdosResultsPanel(ResultsPanel[PdosResultsModel]):
    def _render(self):
        pdos_node = self._model.fetch_child_process_node()
        model = BandsPdosModel.from_nodes(pdos_node=pdos_node)
        widget = BandsPdosWidget(model=model)
        widget.render()
        self.children = [widget]
