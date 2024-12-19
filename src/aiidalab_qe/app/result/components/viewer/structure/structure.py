from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_widgets_base.viewers import StructureDataViewer

from .model import StructureResultsModel


class StructureResults(ResultsPanel[StructureResultsModel]):
    title = "Final Geometry"
    identifier = "structure"

    def _render(self):
        if not hasattr(self, "widget"):
            self.widget = StructureDataViewer(structure=self._model.outputs.structure)
            self.children = [self.widget]

        # HACK to resize the NGL viewer in cases where it auto-rendered when its
        # container was not displayed, which leads to a null width. This hack restores
        # the original dimensions.
        ngl = self.widget._viewer
        ngl._set_size("100%", "300px")
        ngl.control.zoom(0.0)
