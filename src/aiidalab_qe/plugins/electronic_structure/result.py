"""Electronic structure results view widgets"""

from aiidalab_qe.common.bandpdoswidget import BandPdosWidget
from aiidalab_qe.common.panel import ResultPanel


class Result(ResultPanel):
    title = "Electronic Structure"
    workchain_labels = ["bands", "pdos"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)

    def _update_view(self):
        """Update the view of the widget."""
        #
        try:
            pdos_node = self.node.outputs.pdos
        except AttributeError:
            pdos_node = None

        try:
            if "bands" in self.node.outputs.bands:
                bands_node = self.node.outputs.bands.bands
            elif "bands_projwfc" in self.node.outputs.bands:
                bands_node = self.node.outputs.bands.bands_projwfc
        except AttributeError:
            bands_node = None
        _bands_dos_widget = BandPdosWidget(bands=bands_node, pdos=pdos_node)
        # update the electronic structure tab
        self.children = [_bands_dos_widget]
