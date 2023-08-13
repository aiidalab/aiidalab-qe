"""Bands results view widgets

"""


from aiidalab_qe.app.common.panel import ResultPanel


def export_bands_data(outputs, fermi_energy=None):
    """Export the bands data from the outputs of the calculation."""
    import json

    from monty.json import jsanitize

    if "band_structure" in outputs:
        data = json.loads(
            outputs.band_structure._exportcontent("json", comments=False)[0]
        )
        # The fermi energy from band calculation is not robust.
        data["fermi_level"] = fermi_energy or outputs.band_parameters["fermi_energy"]
        return [
            jsanitize(data),
        ]
    else:
        return None


class Result(ResultPanel):
    """Result panel for the bands calculation."""

    title = "Bands"

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, identifier="bands", **kwargs)

    def _update_view(self):
        from widget_bandsplot import BandsPlotWidget

        bands_data = export_bands_data(self.outputs)
        _bands_plot_view = BandsPlotWidget(
            bands=bands_data,
        )
        self.children = [
            _bands_plot_view,
        ]
