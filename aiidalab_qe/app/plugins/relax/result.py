"""Relax results view widgets

"""

from aiidalab_qe.app.panel import ResultPanel


class Result(ResultPanel):
    title = "Final Geometry"
    workchain_label = "relax"

    def _update_view(self):
        from aiidalab_widgets_base.viewers import StructureDataViewer

        self._structure_view = StructureDataViewer(
            structure=self.node.outputs.output_structure
        )
        self.children = [self._structure_view]
