from aiidalab_qe.common.panel import ResultsModel


class StructureResultsModel(ResultsModel):
    identifier = "structure"

    include = True
    _this_process_label = "PwRelaxWorkChain"
