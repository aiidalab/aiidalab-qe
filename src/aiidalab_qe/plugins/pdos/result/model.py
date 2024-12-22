from aiidalab_qe.common.panel import ResultsModel


class PdosResultsModel(ResultsModel):
    title = "PDOS"
    identifier = "pdos"

    _this_process_label = "PdosWorkChain"

    def get_pdos_node(self):
        return self._get_child_outputs()
