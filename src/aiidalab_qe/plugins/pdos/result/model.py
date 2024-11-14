from aiida.common.extendeddicts import AttributeDict
from aiidalab_qe.common.panel import ResultsModel


class PdosResultsModel(ResultsModel):
    identifier = "pdos"

    _this_process_label = "PdosWorkChain"

    def get_pdos_node(self):
        if not (node := self._fetch_child_process_node()):
            return
        outputs = {key: getattr(node.outputs, key) for key in node.outputs}
        return AttributeDict(outputs)

    @property
    def _has_dos(self):
        node = self._fetch_child_process_node()
        return node and "dos" in node.outputs

    @property
    def _has_dos_projections(self):
        node = self._fetch_child_process_node()
        return node and "projwfc" in node.outputs
