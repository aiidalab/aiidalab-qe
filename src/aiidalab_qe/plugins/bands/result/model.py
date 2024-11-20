from aiidalab_qe.common.panel import ResultsModel


class BandsResultsModel(ResultsModel):
    identifier = "bands"

    _this_process_label = "BandsWorkChain"

    def get_bands_node(self):
        outputs = self._get_child_outputs()
        if "bands" in outputs:
            return outputs.bands
        elif "bands_projwfc" in outputs:
            return outputs.bands_projwfc
        else:
            # If neither 'bands' nor 'bands_projwfc' exist, use 'bands_output' itself
            # This is the case for compatibility with older versions of the plugin
            return outputs
