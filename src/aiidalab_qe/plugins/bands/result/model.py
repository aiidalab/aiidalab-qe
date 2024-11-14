from aiidalab_qe.common.panel import ResultsModel


class BandsResultsModel(ResultsModel):
    identifier = "bands"

    _this_process_label = "BandsWorkChain"

    def get_bands_node(self):
        if not (node := self._fetch_child_process_node()):
            return None
        if self._has_bands:
            return node.outputs.bands
        elif self._has_band_projections:
            return node.outputs.bands_projwfc
        else:
            # If neither 'bands' nor 'bands_projwfc' exist, use 'bands_output' itself
            # This is the case for compatibility with older versions of the plugin
            return node.outputs

    @property
    def _has_bands(self):
        node = self._fetch_child_process_node()
        return node and "bands" in node.outputs

    @property
    def _has_band_projections(self):
        node = self._fetch_child_process_node()
        return node and "bands_projwfc" in node.outputs
