"""Relax results view widgets

"""

from aiidalab_qe.app.common.panel import ResultPanel


class Result(ResultPanel):
    title = "Final Geometry"

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, identifier="relax", **kwargs)

    def _update_view(self):
        """Update the view with the output structure."""
        from aiidalab_widgets_base.viewers import StructureDataViewer

        self._structure_view = StructureDataViewer(
            structure=self.outputs.output_structure
        )
        self.children = [self._structure_view]
