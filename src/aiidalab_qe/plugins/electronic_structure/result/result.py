"""Electronic structure results view widgets"""

from aiidalab_qe.common.bands_pdos import BandsPdosModel, BandsPdosWidget
from aiidalab_qe.common.panel import ResultsPanel

from .model import ElectronicStructureResultsModel


class ElectronicStructureResultsPanel(ResultsPanel[ElectronicStructureResultsModel]):
    def _render(self):
        bands_node = self._model.fetch_child_process_node("bands")
        pdos_node = self._model.fetch_child_process_node("pdos")
        model = BandsPdosModel.from_nodes(bands_node=bands_node, pdos_node=pdos_node)
        widget = BandsPdosWidget(model=model)
        widget.render()
        self.children = [widget]
