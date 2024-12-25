"""PDOS results view widgets"""

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import PdosResultsModel


class PdosResultsPanel(ResultsPanel[PdosResultsModel]):
    def _render(self):
        pdos_node = self._model.get_pdos_node()
        model = BandsPdosModel()
        widget = BandsPdosWidget(model=model, pdos=pdos_node)
        widget.render()
        self.children = [widget]
