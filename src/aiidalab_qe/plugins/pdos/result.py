"""PDOS results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.bandpdoswidget import BandPdosWidget
from aiidalab_qe.common.panel import ResultPanel


class Result(ResultPanel):
    title = "PDOS"
    workchain_labels = ["pdos"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)

    def _update_view(self):
        """Update the view of the widget."""

        try:
            pdos_node = self.node.outputs.pdos
        except AttributeError:
            pdos_node = None

        _pdos_plot_view = BandPdosWidget(pdos=pdos_node)

        _pdos_plot_view.layout = ipw.Layout(width="1000px")

        self.children = [_pdos_plot_view]
