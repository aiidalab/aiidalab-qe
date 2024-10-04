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

        try:
            bands_node = self.node.outputs.bands
        except AttributeError:
            bands_node = None
        _bands_dos_widget = BandPdosWidget(bands=bands_node, pdos=pdos_node)

        plot_container = ipw.HTML(
            """
            <div style='overflow-x: auto; width: 80%; max-width: 80%;'>
                <!-- Placeholder for the Plotly plot widget -->
            </div>
            """
        )
        _bands_dos_widget.layout = ipw.Layout(width="1000px")
        self.children = [
            ipw.VBox(
                [
                    plot_container,  # The scrollable container
                    _bands_dos_widget,  # The actual Plotly widget
                ]
            )
        ]
