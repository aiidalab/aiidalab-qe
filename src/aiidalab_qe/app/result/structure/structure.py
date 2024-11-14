from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_widgets_base.viewers import StructureDataViewer

from .model import StructureResultsModel


class StructureResults(ResultsPanel[StructureResultsModel]):
    title = "Final Geometry"
    identifier = "structure"

    def render(self):
        if self.rendered:
            return
        widget = StructureDataViewer(structure=self._model.outputs.structure)
        self.children = [widget]
        self.rendered = True
