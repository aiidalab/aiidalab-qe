"""Bands results view widgets"""

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import BandsResultsModel


class BandsResultsPanel(ResultsPanel[BandsResultsModel]):
    def _render(self):
        bands_node = self._model.fetch_child_process_node()
        model = BandsPdosModel.from_nodes(bands=bands_node)
        widget = BandsPdosWidget(model=model)
        widget.render()
        self.children = [widget]
