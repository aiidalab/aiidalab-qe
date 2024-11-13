from aiidalab_qe.common.panel import ResultsModel


class PdosResultsModel(ResultsModel):
    identifier = "pdos"

    _this_process_label = "PdosWorkChain"

    def get_pdos_node(self):
        try:
            return self.outputs.pdos
        except AttributeError:
            return None

    @property
    def _has_dos(self):
        return self.has_results and hasattr(self.outputs.pdos, "dos")

    @property
    def _has_dos_projections(self):
        return self.has_results and hasattr(self.outputs.pdos, "projwfc")
