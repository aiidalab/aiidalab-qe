"""Electronic structure results view widgets"""

import ipywidgets as ipw

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

        _bands_dos_widget = BandPdosWidget(bands=bands_node, pdos=pdos_node)

        _bands_dos_widget.layout = ipw.Layout(width="1000px")
        self.children = [_bands_dos_widget]
