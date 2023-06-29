"""Bands results view widgets

"""


from aiidalab_qe.app.panel import ResultPanel


def export_bands_data(node, fermi_energy=None):
    import json

    from monty.json import jsanitize

    if "band_structure" in node.outputs:
        data = json.loads(
            node.outputs.band_structure._exportcontent("json", comments=False)[0]
        )
        # The fermi energy from band calculation is not robust.
        data["fermi_level"] = (
            fermi_energy or node.outputs.band_parameters["fermi_energy"]
        )
        return [
            jsanitize(data),
        ]
    else:
        return None


class Result(ResultPanel):
    title = "Bands"

    def _update_view(self):
        from widget_bandsplot import BandsPlotWidget

        bands_data = export_bands_data(self.node)
        _bands_plot_view = BandsPlotWidget(
            bands=bands_data,
        )
        self.children = [
            _bands_plot_view,
        ]
