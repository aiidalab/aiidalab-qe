"""Bands results view widgets

"""

from aiidalab_widgets_base import register_viewer_widget

from aiidalab_qe.common.panel import ResultPanel


def export_bands_data(work_chain_node, fermi_energy=None):
    """Export the bands data from the outputs of the calculation."""
    import json

    from monty.json import jsanitize

    if "band_structure" in work_chain_node.outputs:
        data = json.loads(
            work_chain_node.outputs.band_structure._exportcontent(
                "json", comments=False
            )[0]
        )
        # The fermi energy from band calculation is not robust.
        data["fermi_level"] = (
            fermi_energy or work_chain_node.outputs.band_parameters["fermi_energy"]
        )
        return [
            jsanitize(data),
        ]
    else:
        return None


@register_viewer_widget("aiida.workflows:quantumespresso.pw.bands")
class Result(ResultPanel):
    """Result panel for the bands calculation."""

    title = "Bands"
    workchain_labels = ["bands"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)
        self.workchain_nodes = {}
        if node.process_type == "aiida.workflows:quantumespresso.pw.bands":
            self.workchain_nodes["bands"] = node
        elif node.process_type == "aiidalab_qe.workflows.QeAppWorkChain":
            if "bands" in node.base.links.get_outgoing().all_link_labels():
                self.workchain_nodes[
                    "bands"
                ] = node.base.links.get_outgoing().get_node_by_label("bands")
            else:
                self.workchain_nodes["bands"] = None
        self._update_view()

    def _update_view(self):
        from widget_bandsplot import BandsPlotWidget

        if self.workchain_nodes["bands"] is None:
            return
        bands_data = export_bands_data(self.workchain_nodes["bands"])
        _bands_plot_view = BandsPlotWidget(
            bands=bands_data,
        )
        self.children = [
            _bands_plot_view,
        ]
