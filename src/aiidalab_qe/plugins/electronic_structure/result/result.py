"""Electronic structure results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.bands_pdos import BandPdosWidget, BandsPdosModel
from aiidalab_qe.common.panel import ResultsPanel

from .model import ElectronicStructureResultModel


class ElectronicStructureResult(ResultsPanel[ElectronicStructureResultModel]):
    title = "Electronic Structure"
    identifier = "electronic_structure"
    workchain_labels = ["bands", "pdos"]

    def render(self):
        if self.rendered:
            return

        pdos_node = self._model.get_pdos_node()
        bands_node = self._model.get_bands_node()

        model = BandsPdosModel()
        widget = BandPdosWidget(model=model, bands=bands_node, pdos=pdos_node)
        widget.layout = ipw.Layout(width="1000px")
        widget.render()

        self.children = [widget]

        self.rendered = True
