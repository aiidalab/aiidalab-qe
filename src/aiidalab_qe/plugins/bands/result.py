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
        # Initialize bands_node to None by default
        bands_node = None

        # Check if the workchain has the 'bands' output
        if hasattr(self.node.outputs, "bands"):
            bands_output = self.node.outputs.bands

            # Check for 'bands' or 'bands_projwfc' attributes within 'bands' output
            if hasattr(bands_output, "bands"):
                bands_node = bands_output.bands
            elif hasattr(bands_output, "bands_projwfc"):
                bands_node = bands_output.bands_projwfc
            else:
                # If neither 'bands' nor 'bands_projwfc' exist, use 'bands_output' itself
                # This is the case for compatibility with older versions of the plugin
                bands_node = bands_output

        _bands_plot_view = BandPdosWidget(bands=bands_node)
        self.children = [
            _bands_plot_view,
        ]
