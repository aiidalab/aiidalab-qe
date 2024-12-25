"""Electronic structure results view widgets"""

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import ElectronicStructureResultsModel


class ElectronicStructureResultsPanel(ResultsPanel[ElectronicStructureResultsModel]):
    def _render(self):
        bands_node = self._model.get_bands_node()
        pdos_node = self._model.get_pdos_node()
        model = BandsPdosModel()
        widget = BandsPdosWidget(model=model, bands=bands_node, pdos=pdos_node)
        widget.render()
        self.children = [widget]
