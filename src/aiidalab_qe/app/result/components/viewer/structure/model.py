from aiidalab_qe.common.panel import ResultsModel


class StructureResultsModel(ResultsModel):
    title = "Structure"
    identifier = "structure"

    _this_process_label = "PwRelaxWorkChain"

    source = None

    @property
    def include(self):
        return True

    def update(self):
        is_relaxed = "relax" in self.properties
        self.title = "Relaxed structure" if is_relaxed else "Initial structure"
        self.source = self.outputs if is_relaxed else self.inputs
        self.auto_render = not is_relaxed or self.has_results

    def get_structure(self):
        try:
            return self.source.structure if self.source else None
        except AttributeError:
            # If source is outputs but job failed, there may not be a structure
            return None
