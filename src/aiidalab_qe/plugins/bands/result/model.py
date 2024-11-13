from aiidalab_qe.common.panel import ResultsModel


class BandsResultsModel(ResultsModel):
    identifier = "bands"

    _this_process_label = "BandsWorkChain"

    def get_bands_node(self):
        # Check for 'bands' or 'bands_projwfc' attributes within 'bands' output
        if self._has_bands:
            return self.outputs.bands.bands
        elif self._has_band_projections:
            return self.outputs.bands.bands_projwfc
        elif self.has_results:
            # If neither 'bands' nor 'bands_projwfc' exist, use 'bands_output' itself
            # This is the case for compatibility with older versions of the plugin
            return self.outputs.bands

    @property
    def _has_bands(self):
        return self.has_results and hasattr(self.outputs.bands, "bands")

    @property
    def _has_band_projections(self):
        return self.has_results and hasattr(self.outputs.bands, "bands_projwfc")
