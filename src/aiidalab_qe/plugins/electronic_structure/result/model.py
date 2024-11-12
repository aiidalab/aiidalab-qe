from aiidalab_qe.common.panel import ResultsModel


class ElectronicStructureResultModel(ResultsModel):
    def get_pdos_node(self):
        try:
            return self.outputs.pdos
        except AttributeError:
            return None

    def get_bands_node(self):
        # Check for 'bands' or 'bands_projwfc' attributes within 'bands' output
        if self._has_bands:
            return self.outputs.bands.bands
        elif self._has_bands_projwfc:
            return self.outputs.bands.bands_projwfc
        elif self._has_bands_output:
            # If neither 'bands' nor 'bands_projwfc' exist, use 'bands_output' itself
            # This is the case for compatibility with older versions of the plugin
            return self.outputs.bands

    @property
    def _has_bands_output(self):
        return hasattr(self.outputs, "bands")

    @property
    def _has_bands(self):
        if not self._has_bands_output:
            return False
        return hasattr(self.outputs.bands, "bands")

    @property
    def _has_bands_projwfc(self):
        if not self._has_bands_output:
            return False
        return hasattr(self.outputs.bands, "bands_projwfc")
