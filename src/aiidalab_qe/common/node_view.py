"""Results view widgets (MOVE TO OTHER MODULE!)

Authors: AiiDAlab team
"""

import ipywidgets as ipw
import nglview
import traitlets as tl
from ase import Atoms

from aiida import orm
from aiidalab_widgets_base import register_viewer_widget

from .widgets import CalcJobOutputFollower, LogOutputWidget


class MinimalStructureViewer(ipw.VBox):
    structure = tl.Union([tl.Instance(Atoms), tl.Instance(orm.Node)], allow_none=True)
    _displayed_structure = tl.Instance(Atoms, allow_none=True, read_only=True)

    background = tl.Unicode()
    supercell = tl.List(tl.Int())

    def __init__(self, structure, *args, **kwargs):
        self._viewer = nglview.NGLWidget()
        self._viewer.camera = "orthographic"
        self._viewer.stage.set_parameters(mouse_preset="pymol")
        ipw.link((self, "background"), (self._viewer, "background"))

        self.structure = structure

        super().__init__(
            *args,
            children=[
                self._viewer,
            ],
            **kwargs,
        )

    @tl.default("background")
    def _default_background(self):
        return "#FFFFFF"

    @tl.default("supercell")
    def _default_supercell(self):
        return [1, 1, 1]

    @tl.validate("structure")
    def _valid_structure(self, change):
        """Update structure."""
        structure = change["value"]

        if structure is None:
            return None  # if no structure provided, the rest of the code can be skipped

        if isinstance(structure, Atoms):
            return structure
        if isinstance(structure, orm.Node):
            return structure.get_ase()
        raise ValueError(
            "Unsupported type {}, structure must be one of the following types: "
            "ASE Atoms object, AiiDA CifData or StructureData."
        )

    @tl.observe("structure")
    def _update_displayed_structure(self, change):
        """Update displayed_structure trait after the structure trait has been modified."""
        # Remove the current structure(s) from the viewer.
        if change["new"] is not None:
            self.set_trait("_displayed_structure", change["new"].repeat(self.supercell))
        else:
            self.set_trait("_displayed_structure", None)

    @tl.observe("_displayed_structure")
    def _update_structure_viewer(self, change):
        """Update the view if displayed_structure trait was modified."""
        with self.hold_trait_notifications():
            for comp_id in self._viewer._ngl_component_ids:
                self._viewer.remove_component(comp_id)
            self.selection = []
            if change["new"] is not None:
                self._viewer.add_component(nglview.ASEStructure(change["new"]))
                self._viewer.clear()
                self._viewer.stage.set_parameters(clipDist=0)
                self._viewer.add_representation("unitcell", diffuse="#df0587")
                self._viewer.add_representation("ball+stick", aspectRatio=3.5)


class VBoxWithCaption(ipw.VBox):
    def __init__(self, caption, body, *args, **kwargs):
        super().__init__(*args, children=[ipw.HTML(caption), body], **kwargs)


@register_viewer_widget("process.calculation.calcjob.CalcJobNode.")
class CalcJobNodeViewerWidget(ipw.VBox):
    def __init__(self, calcjob, **kwargs):
        self.calcjob = calcjob
        self.output_follower = CalcJobOutputFollower()
        self.log_output = LogOutputWidget()

        self.output_follower.calcjob_uuid = self.calcjob.uuid
        self.output_follower.observe(self._observe_output_follower_lineno, ["lineno"])

        super().__init__(
            children=[
                ipw.HTML(f"CalcJob: {self.calcjob}"),
                self.log_output,
            ],
            layout=ipw.Layout(height="100%"),
            **kwargs,
        )
        self.add_class("calcjob-node-viewer")

    def _observe_output_follower_lineno(self, _):
        with self.hold_trait_notifications():
            self.log_output.filename = self.output_follower.filename
            self.log_output.value = "\n".join(self.output_follower.output)
