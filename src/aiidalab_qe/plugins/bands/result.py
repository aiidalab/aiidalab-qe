"""Bands results view widgets"""

from aiidalab_qe.common.bandpdoswidget import BandPdosWidget
from aiidalab_qe.common.panel import ResultPanel


class Result(ResultPanel):
    """Result panel for the bands calculation."""

    title = "Bands"
    workchain_labels = ["bands"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)

    def _update_view(self):
        # Check if the workchain has the outputs
        try:
            if "bands" in self.node.outputs.bands:
                bands_node = self.node.outputs.bands.bands
            elif "bands_projwfc" in self.node.outputs.bands:
                bands_node = self.node.outputs.bands.bands_projwfc
        except AttributeError:
            bands_node = None

        _bands_plot_view = BandPdosWidget(bands=bands_node)
        self.children = [
            _bands_plot_view,
        ]
