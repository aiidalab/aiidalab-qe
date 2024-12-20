from aiidalab_qe.common.panel import ResultsModel


class StructureResultsModel(ResultsModel):
    identifier = "structure"

    _this_process_label = "PwRelaxWorkChain"

    @property
    def include(self):
        return "relax" in self.properties
