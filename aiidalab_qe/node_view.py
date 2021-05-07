"""Results view widgets (MOVE TO OTHER MODULE!)

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""

import json

import ipywidgets as ipw
import nglview
from aiida.common import NotExistentAttributeError
from aiida.orm import Node
from aiidalab_widgets_base import register_viewer_widget
from aiidalab_widgets_base import ProcessMonitor
from ase import Atoms
from jinja2 import Environment
from monty.json import MontyEncoder
from monty.json import jsanitize
from traitlets import Instance
from traitlets import Int
from traitlets import List
from traitlets import Unicode
from traitlets import Union
from traitlets import default
from traitlets import observe
from traitlets import validate
from widget_bandsplot import BandsPlotWidget

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
            return "yes" if truthy else "no"

        style = "background-color: #fafafa; line-height: normal"
        report = generate_report_dict(self.wc_node)

        env = Environment()
        env.filters.update(
            {
                "fmt_yes_no": _fmt_yes_no,
            }
        )
        template = env.from_string(
            """
            <pre style="{style}">
            <table>
                <tr>
                    <td>Structure geometry optimized:</td>
                    <td>{{ relaxed | fmt_yes_no }}</td>
                </tr>
                <tr>
                    <td>Pseudopotential library:</td>
                    <td>{{ pseudo_family }}</td>
                </tr>
                <tr>
                    <td>Energy cutoff (wave functions):</td>
                    <td>{{ energy_cutoff_wfc }}</td>
                </tr>
                <tr>
                    <td>Energy cutoff (charge density):</td>
                    <td>{{ energy_cutoff_rho }}</td>
                </tr>
                <tr>
                    <td>K-point mesh distance (SCF):</td>
                    <td>{{ scf_kpoints_distance }}</td>
                </tr>
                {% if bands_computed %}
                <tr>
                    <td>K-point line distance (Bands)</td>
                    <td>{{ bands_kpoints_distance }}</td>
                </tr>
                {% endif %}
                {% if pdos_computed %}
                <tr>
                    <td>K-point mesh distance (NSCF)</td>
                    <td>{{ nscf_kpoints_distance }}</td>
                </tr>
                {% endif %}
            </table>
            </pre>
            """
        )

        self.summary_view = ipw.HTML(template.render(style=style, **report))

        super().__init__(children=[self.summary_view], **kwargs)


@register_viewer_widget("process.workflow.workchain.WorkChainNode.")
class WorkChainViewer(ipw.VBox):
    def __init__(self, node, **kwargs):
        if node.process_label != "QeAppWorkChain":
            raise KeyError(str(node.node_type))

        self.node = node

        self.output_structure_view = None
        self.output_structure_view_box = ipw.HBox(
            children=[ipw.HTML("No structure data available yet.")],
            layout=ipw.Layout(
                flex="1 1 auto",
                min_width="200px",
                min_height="200px",
                max_width="300px",
                max_height="300px",
                border="1px solid grey",
            ),
        )

        self.bands_plot = None
        self.bands_plot_box = ipw.HBox(
            children=[ipw.HTML("No bands data available yet.")],
            layout=ipw.Layout(
                min_height="300px",
                flex="2 1 auto",
            ),
        )

        structure_formula = self.node.inputs.structure.get_formula()
        self.title = ipw.HTML(
            f"<h4>QE App Work Chain (pk: {self.node.pk}) &mdash; {structure_formula} </h4>"
        )
        self.summary_view = SummaryView(
            self.node,
            layout=ipw.Layout(
                flex="1 1 auto",
            ),
        )

        self._process_monitor = ProcessMonitor(
            process=self.node, callbacks=[self._update_view]
        )

        super().__init__(
            children=[
                self.title,
                ipw.HBox(
                    [
                        ipw.VBox(
                            [self.output_structure_view_box, self.summary_view],
                            layout=ipw.Layout(flex="1 1 auto"),
                        ),
                        self.bands_plot_box,
                    ],
                ),
            ],
            **kwargs,
        )

    def _update_view(self):
        if self.output_structure_view is None:
            try:
                self.output_structure_view = MinimalStructureViewer(
                    structure=self.node.outputs.structure,
                    layout=ipw.Layout(flex="1 1 auto"),
                )
                self.output_structure_view_box.children = [self.output_structure_view]
            except NotExistentAttributeError:
                pass

        if self.bands_plot is None:
            bands_data = export_bands_data(self.node)
            pdos_data = export_pdos_data(self.node)

            if bands_data or pdos_data:
                self.bands_plot = BandsPlotWidget(
                    bands=[bands_data] if bands_data else [],
                    dos=pdos_data,
                    plot_fermilevel=True,
                )
                self.bands_plot_box.children = [self.bands_plot]
