"""Results view widgets (MOVE TO OTHER MODULE!)

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""

import json

import ipywidgets as ipw
from aiidalab_widgets_base import register_viewer_widget
from aiidalab_widgets_base import viewer
from monty.json import MontyEncoder
from monty.json import jsanitize
from widget_bandsplot import BandsPlotWidget


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
    def __init__(self, inputs, **kwargs):
        self.inputs = inputs

        def _fmt_truthy(truthy):
            return "yes" if truthy else "false"

        style = "background-color: #fafafa; line-height: normal"

        self.summary_view = ipw.HTML(
            f"""
            <pre style="{style}">
            <table>
                <tr>
                    <td>Structure relaxation:</td>
                    <td>{_fmt_truthy('relax__base__pw__parameters' in inputs)}</td>
                </tr>
                <tr>
                    <td>Bands computed:</td>
                    <td>{_fmt_truthy('bands__bands__pw__parameters' in inputs)}</td>
                </tr>
                <tr>
                    <td>Density of State (DoS) computed:</td>
                    <td>{_fmt_truthy('pdos__dos__parameters' in inputs)}</td>
                </tr>
            </table>
            </pre>
            """
        )

        super().__init__(children=[self.summary_view], **kwargs)


@register_viewer_widget("process.workflow.workchain.WorkChainNode.")
class WorkChainViewer(ipw.VBox):
    def __init__(self, node, **kwargs):
        self.node = node

        self.output_structure_view = viewer(
            self.node.outputs.structure,
            configure_view=False,
            layout=ipw.Layout(flex="1 0 auto"),
        )

        bands_data = export_bands_data(self.node)
        pdos_data = export_pdos_data(self.node)

        if bands_data or pdos_data:
            self.bands_plot = BandsPlotWidget(
                bands=[bands_data] if bands_data else [],
                dos=pdos_data,
                plot_fermilevel=True,
            )
        else:
            self.bands_plot = ipw.HTML("No bands or density of state data computed.")

        bands_plot_box = ipw.HBox(
            children=[self.bands_plot],
            layout=ipw.Layout(
                min_height="300px",
            ),
        )

        self.title = ipw.HTML(f"<h4>QE App Work Chain (pk: {self.node.pk})</h4>")
        self.summary_view = SummaryView(node.inputs, layout=ipw.Layout(flex="1 0 auto"))

        super().__init__(
            children=[
                self.title,
                ipw.HBox(
                    [
                        self.summary_view,
                        self.output_structure_view,
                    ]
                ),
                bands_plot_box,
            ],
            **kwargs,
        )
