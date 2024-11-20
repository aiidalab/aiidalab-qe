"""Bands results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import BandsResultsModel


class BandsResults(ResultsPanel[BandsResultsModel]):
    title = "Bands"
    identifier = "bands"
    workchain_labels = ["bands"]

    def render(self):
        if self.rendered:
            return
        bands_node = self._model.get_bands_node()
        model = BandsPdosModel()
        widget = BandsPdosWidget(model=model, bands=bands_node)
        widget.layout = ipw.Layout(width="1000px")
        widget.render()
        self.children = [widget]
        self.rendered = True
