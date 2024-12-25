"""Bands results view widgets"""

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import BandsResultsModel


class BandsResultsPanel(ResultsPanel[BandsResultsModel]):
    def _render(self):
        bands_node = self._model.get_bands_node()
        model = BandsPdosModel()
        widget = BandsPdosWidget(model=model, bands=bands_node)
        widget.render()
        self.children = [widget]
