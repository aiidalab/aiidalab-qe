"""Results view widgets (MOVE TO OTHER MODULE!)

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""

import json
from importlib import resources

import ipywidgets as ipw
import nglview
from aiida.orm import Node
from aiidalab_widgets_base import ProcessMonitor, register_viewer_widget
from aiidalab_widgets_base.viewers import StructureDataViewer
from ase import Atoms
from jinja2 import Environment
from monty.json import MontyEncoder, jsanitize
from traitlets import Instance, Int, List, Unicode, Union, default, observe, validate
from widget_bandsplot import BandsPlotWidget

from aiidalab_qe import static
from aiidalab_qe.report import generate_report_dict


class MinimalStructureViewer(ipw.VBox):

    structure = Union([Instance(Atoms), Instance(Node)], allow_none=True)
    _displayed_structure = Instance(Atoms, allow_none=True, read_only=True)

    background = Unicode()
    supercell = List(Int)

    def __init__(self, structure, *args, **kwargs):

        self._viewer = nglview.NGLWidget()
        self._viewer.camera = "orthographic"
        self._viewer.stage.set_parameters(mouse_preset="pymol")
        ipw.link((self, "background"), (self._viewer, "background"))

        self.structure = structure

        super().__init__(
            children=[
                self._viewer,
            ],
            *args,
            **kwargs,
        )

    @default("background")
    def _default_background(self):
        return "#FFFFFF"

    @default("supercell")
    def _default_supercell(self):
        return [1, 1, 1]

    @validate("structure")
    def _valid_structure(self, change):  # pylint: disable=no-self-use
        """Update structure."""
        structure = change["value"]

        if structure is None:
            return None  # if no structure provided, the rest of the code can be skipped

        if isinstance(structure, Atoms):
            return structure
        if isinstance(structure, Node):
            return structure.get_ase()
        raise ValueError(
            "Unsupported type {}, structure must be one of the following types: "
            "ASE Atoms object, AiiDA CifData or StructureData."
        )

    @observe("structure")
    def _update_displayed_structure(self, change):
        """Update displayed_structure trait after the structure trait has been modified."""
        # Remove the current structure(s) from the viewer.
        if change["new"] is not None:
            self.set_trait("_displayed_structure", change["new"].repeat(self.supercell))
        else:
            self.set_trait("_displayed_structure", None)

    @observe("_displayed_structure")
    def _update_structure_viewer(self, change):
        """Update the view if displayed_structure trait was modified."""
        with self.hold_trait_notifications():
            for (
                comp_id
            ) in self._viewer._ngl_component_ids:  # pylint: disable=protected-access
                self._viewer.remove_component(comp_id)
            self.selection = list()
            if change["new"] is not None:
                self._viewer.add_component(nglview.ASEStructure(change["new"]))
                self._viewer.clear()
                self._viewer.stage.set_parameters(clipDist=0)
                self._viewer.add_representation("unitcell", diffuse="#df0587")
                self._viewer.add_representation("ball+stick", aspectRatio=3.5)


def export_bands_data(work_chain_node):
    if "band_structure" in work_chain_node.outputs:
        data = json.loads(
            work_chain_node.outputs.band_structure._exportcontent(
                "json", comments=False
            )[0]
        )
        data["fermi_level"] = work_chain_node.outputs.band_parameters["fermi_energy"]
        return jsanitize(data)


def export_pdos_data(work_chain_node):
    if "dos" in work_chain_node.outputs:
        fermi_energy = work_chain_node.outputs.band_parameters["fermi_energy"]
        x_label, energy_dos, energy_units = work_chain_node.outputs.dos.get_x()
        tdos_values = {
            f"{n} | {u}": v for n, v, u in work_chain_node.outputs.dos.get_y()
        }

        pdos_orbitals = []

        for orbital, pdos, energy in work_chain_node.outputs.projections.get_pdos():
            orbital_data = orbital.get_orbital_dict()
            kind_name = orbital_data["kind_name"]
            orbital_name = orbital.get_name_from_quantum_numbers(
                orbital_data["angular_momentum"], orbital_data["magnetic_number"]
            )
            pdos_orbitals.append(
                {
                    "kind": kind_name,
                    "orbital": orbital_name,
                    "energy | eV": energy,
                    "pdos | states/eV": pdos,
                }
            )

        data_dict = {
            "fermi_energy": fermi_energy,
            "tdos": {f"energy | {energy_units}": energy_dos, "values": tdos_values},
            "pdos": pdos_orbitals,
        }

        # And this is why we shouldn't use special encoders...
        return json.loads(json.dumps(data_dict, cls=MontyEncoder))

        fermi_energy = work_chain_node.outputs.nscf__output_parameters.get_dict()[
            "fermi_energy"
        ]
        (
            xlabel,
            energy_dos,
            energy_units,
        ) = work_chain_node.outputs.dos__output_dos.get_x()
        tdos_values = {
            f"{n} | {u}": v
            for n, v, u in work_chain_node.outputs.dos__output_dos.get_y()
        }


def export_data(work_chain_node):
    return dict(
        bands=export_bands_data(work_chain_node), dos=export_pdos_data(work_chain_node)
    )


class VBoxWithCaption(ipw.VBox):
    def __init__(self, caption, body, *args, **kwargs):
        super().__init__(children=[ipw.HTML(caption), body], *args, **kwargs)


class SummaryView(ipw.VBox):
    def __init__(self, wc_node, **kwargs):

        self.wc_node = wc_node

        def _fmt_yes_no(truthy):
            return "Yes" if truthy else "No"

        report = generate_report_dict(self.wc_node)

        env = Environment()
        env.filters.update(
            {
                "fmt_yes_no": _fmt_yes_no,
            }
        )
        template = resources.read_text(static, "workflow_summary.jinja")
        style = resources.read_text(static, "style.css")
        self.summary_view = ipw.HTML(
            env.from_string(template).render(style=style, **report)
        )
        super().__init__(
            children=[self.summary_view],
            layout=ipw.Layout(min_height="380px"),
            **kwargs,
        )


@register_viewer_widget("process.workflow.workchain.WorkChainNode.")
class WorkChainViewer(ipw.VBox):
    def __init__(self, node, **kwargs):
        if node.process_label != "QeAppWorkChain":
            raise KeyError(str(node.node_type))

        self.node = node
        self.status = []

        structure_formula = self.node.inputs.structure.get_formula()
        self.title = ipw.HTML(
            f"""
            <hr style="height:2px;background-color:#2097F3;">
            <h4>QE App Work Chain (pk: {self.node.pk}) &mdash; {structure_formula} </h4>
        """
        )

        self.result_titles = [
            "Workflow Summary",
        ]
        self.result_children = [self.get_summary()]
        self.result_tabs = ipw.Tab()

        # An ugly fix to the structure appearance problem
        # https://github.com/aiidalab/aiidalab-qe/issues/69
        def on_change(change):
            index = change["new"]
            # Accessing the viewer only if the corresponding tab is present.
            if self.result_tabs._titles[str(index)] == "Final Geometry":
                self.struct_view._viewer.handle_resize()

                def toggle_camera():
                    """Toggle camera between perspective and orthographic."""
                    self.struct_view._viewer.camera = (
                        "perspective"
                        if self.struct_view._viewer.camera == "orthographic"
                        else "orthographic"
                    )

                toggle_camera()
                toggle_camera()

        self.result_tabs.observe(on_change, "selected_index")
        self._update_view(first_run=True)

        super().__init__(
            children=[
                self.title,
                self.result_tabs,
            ],
            **kwargs,
        )
        self._process_monitor = ProcessMonitor(
            process=self.node,
            callbacks=[
                self._update_view,
            ],
        )

    def _update_view(self, first_run=False):

        update = False

        if "structure" not in self.status:
            if "structure" in self.node.outputs:
                self.result_titles.append("Final Geometry")
                self.struct_view = self.get_structure()
                self.result_children.append(self.struct_view)
                self.status.append("structure")
                update = True

        if "bands" not in self.status:
            if "band_structure" in self.node.outputs:
                self.result_titles.append("Band Structure")
                self.result_children.append(self.get_bands())
                self.status.append("bands")
                update = True

        if first_run or update:
            index = self.result_tabs.selected_index

            self.result_tabs.children = self.result_children
            for tab_number, title in enumerate(self.result_titles):
                self.result_tabs.set_title(tab_number, title)

            self.result_tabs.selected_index = index

    def get_summary(self):
        return SummaryView(self.node)

    def get_structure(self):
        return StructureDataViewer(structure=self.node.outputs.structure)

    def get_bands(self):
        bands_data = export_bands_data(self.node)
        if bands_data:
            return BandsPlotWidget(bands=[bands_data], plot_fermilevel=True)
