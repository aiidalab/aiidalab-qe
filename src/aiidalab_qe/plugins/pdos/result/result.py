"""PDOS results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import PdosResultsModel


class PdosResults(ResultsPanel[PdosResultsModel]):
    title = "PDOS"
    identifier = "pdos"
    workchain_labels = ["pdos"]

    def render(self):
        if self.rendered:
            return
        pdos_node = self._model.get_pdos_node()
        model = BandsPdosModel()
        widget = BandsPdosWidget(model=model, pdos=pdos_node)
        widget.layout = ipw.Layout(width="1000px")
        widget.render()
        self.children = [widget]
        self.rendered = True
