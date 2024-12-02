"""Electronic structure results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import ElectronicStructureResultsModel


class ElectronicStructureResultsPanel(ResultsPanel[ElectronicStructureResultsModel]):
    title = "Electronic Structure"
    identifier = "electronic_structure"
    workchain_labels = ["bands", "pdos"]

    def render(self):
        if self.rendered:
            return
        bands_node = self._model.get_bands_node()
        pdos_node = self._model.get_pdos_node()
        model = BandsPdosModel()
        widget = BandsPdosWidget(model=model, bands=bands_node, pdos=pdos_node)
        widget.layout = ipw.Layout(width="1000px")
        widget.render()
        self.children = [widget]
        self.rendered = True
