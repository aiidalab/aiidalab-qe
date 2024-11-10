from aiidalab_qe.common.panel import ResultsModel


class PdosResultModel(ResultsModel):
    def get_pdos_node(self):
        try:
            return self.process_node.outputs.pdos
        except AttributeError:
            return None
