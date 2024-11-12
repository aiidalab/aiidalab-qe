from aiidalab_qe.common.panel import ResultsModel


class PdosResultsModel(ResultsModel):
    def get_pdos_node(self):
        try:
            return self.outputs.pdos
        except AttributeError:
            return None
