import base64
import json
import re

import ipywidgets as ipw
import numpy as np
import plotly.graph_objects as go
from IPython.display import clear_output, display
from plotly.subplots import make_subplots

from aiida.orm import ProjectionData
from aiidalab_widgets_base.utils import StatusHTML, string_range_to_list


class BandPdosPlotly:
    SETTINGS = {
        "axis_linecolor": "#111111",
        "bands_linecolor": "#111111",
        "bands_up_linecolor": "rgba(205, 0, 0, 0.4)",  # Red Opacitiy 40%
        "bands_down_linecolor": "rgba(72,118,255, 0.4)",  # Blue Opacitiy 40%
        "combined_plot_height": 600,
        "combined_plot_width": 900,
        "combined_column_widths": [0.7, 0.3],
        "bands_plot_height": 600,
        "bands_plot_width": 850,
        "pdos_plot_height": 600,
        "pdos_plot_width": 850,
        "vertical_linecolor": "#111111",
        "horizontal_linecolor": "#111111",
        "vertical_range_bands": [-10, 10],
        "horizontal_range_pdos": [-10, 10],
    }

    def __init__(self, bands_data=None, pdos_data=None, project_bands=False):
        self.bands_data = bands_data
        self.pdos_data = pdos_data
        self.fermi_energy = self._get_fermi_energy()
        self.project_bands = project_bands and "projected_bands" in self.bands_data

        # Plotly Axis
        # Plotly settings
        self._bands_xaxis = self._band_xaxis()
        self._bands_yaxis = self._band_yaxis()
        self._pdos_xaxis = self._pdos_xaxis()
        self._pdos_yaxis = self._pdos_yaxis()

    @property
    def plot_type(self):
        """Define the plot type."""
        if self.bands_data and self.pdos_data:
            return "combined"
        elif self.bands_data:
            return "bands"
        elif self.pdos_data:
            return "pdos"

    @property
    def bandspdosfigure(self):
        return self._get_bandspdos_plot()

    def _get_fermi_energy(self):
        """Function to return the Fermi energy information depending on the data available."""
        fermi_data = {}
        if self.pdos_data:
            if "fermi_energy_up" in self.pdos_data:
                fermi_data["fermi_energy_up"] = self.pdos_data["fermi_energy_up"]
                fermi_data["fermi_energy_down"] = self.pdos_data["fermi_energy_down"]
            else:
                fermi_data["fermi_energy"] = self.pdos_data["fermi_energy"]
        else:
            if "fermi_energy_up" in self.bands_data:
                fermi_data["fermi_energy_up"] = self.bands_data["fermi_energy_up"]
                fermi_data["fermi_energy_down"] = self.bands_data["fermi_energy_down"]
            else:
                fermi_data["fermi_energy"] = self.bands_data["fermi_energy"]
        return fermi_data

    def _band_xaxis(self):
        """Function to return the xaxis for the bands plot."""

        if not self.bands_data:
            return None
        paths = self.bands_data.get("paths")
        slider_bands = go.layout.xaxis.Rangeslider(
            thickness=0.08,
            range=[0, paths[-1]["x"][-1]],
        )
        bandxaxis = go.layout.XAxis(
            title="k-points",
            range=[0, paths[-1]["x"][-1]],
            showgrid=True,
            showline=True,
            tickmode="array",
            rangeslider=slider_bands,
            fixedrange=False,
            tickvals=self.bands_data["pathlabels"][1],  # ,self.band_labels[1],
            ticktext=self.bands_data["pathlabels"][0],  # self.band_labels[0],
            showticklabels=True,
            linecolor=self.SETTINGS["axis_linecolor"],
            mirror=True,
            linewidth=2,
            type="linear",
        )

        return bandxaxis

    def _band_yaxis(self):
        """Function to return the yaxis for the bands plot."""

        if not self.bands_data:
            return None

        bandyaxis = go.layout.YAxis(
            title={"text": "Electronic Bands (eV)", "standoff": 1},
            side="left",
            showgrid=True,
            showline=True,
            zeroline=True,
            range=self.SETTINGS["vertical_range_bands"],
            fixedrange=False,
            automargin=True,
            ticks="inside",
            linewidth=2,
            linecolor=self.SETTINGS["axis_linecolor"],
            tickwidth=2,
            zerolinewidth=2,
        )

        return bandyaxis

    def _pdos_xaxis(self):
        """Function to return the xaxis for the pdos plot."""

        if not self.pdos_data:
            return None
        # For combined plot
        axis_settings = {
            "showgrid": True,
            "showline": True,
            "mirror": "ticks",
            "ticks": "inside",
            "linewidth": 2,
            "tickwidth": 2,
            "linecolor": self.SETTINGS["axis_linecolor"],
            "title": "Density of states",
            "side": "bottom",
            "automargin": True,
        }

        if self.plot_type != "combined":
            axis_settings["title"] = "Density of states (eV)"
            axis_settings["range"] = self.SETTINGS["horizontal_range_pdos"]
            axis_settings.pop("side")
            axis_settings.pop("automargin")

        return go.layout.XAxis(**axis_settings)

    def _pdos_yaxis(self):
        """Function to return the yaxis for the pdos plot."""

        if not self.pdos_data:
            return None

        axis_settings = {
            "showgrid": True,
            "showline": True,
            "side": "right" if self.plot_type == "combined" else "left",
            "mirror": "ticks",
            "ticks": "inside",
            "linewidth": 2,
            "tickwidth": 2,
            "linecolor": self.SETTINGS["axis_linecolor"],
            "zerolinewidth": 2,
        }

        return go.layout.YAxis(**axis_settings)

    def _get_bandspdos_plot(self):
        """Function to return the bands plot widget."""

        fig = self._create_fig()
        if self.bands_data:
            self._add_band_traces(fig)

            band_labels = self.bands_data.get("pathlabels")
            for label in band_labels[1]:
                fig.add_vline(
                    x=label,
                    line={"color": self.SETTINGS["vertical_linecolor"], "width": 1},
                )

            if self.project_bands:
                self._add_projection_traces(fig)

        if self.pdos_data:
            self._add_pdos_traces(fig)
            if self.plot_type == "pdos":
                fig.add_vline(
                    x=0,
                    line={
                        "color": self.SETTINGS["vertical_linecolor"],
                        "width": 1,
                        "dash": "dot",
                    },
                )

        if self.plot_type == "combined":
            self._customize_combined_layout(fig)
        else:
            self._customize_single_layout(fig)

        return go.FigureWidget(fig)

    def _create_fig(self):
        """Create a plotly figure.

        The figure layout is different depending on the plot type.
        """
        if self.plot_type != "combined":
            return go.Figure()

        fig = make_subplots(
            rows=1,
            cols=2,
            shared_yaxes=True,
            column_widths=self.SETTINGS["combined_column_widths"],
            horizontal_spacing=0.015,
        )
        return fig

    def _add_traces_to_fig(self, fig, traces, col):
        """Add a list of traces to a figure."""
        if self.plot_type == "combined":
            rows = [1] * len(traces)
            cols = [col] * len(traces)
            fig.add_traces(traces, rows=rows, cols=cols)
        else:
            fig.add_traces(traces)

    def _add_band_traces(self, fig):
        """Generate the band traces and add them to the figure."""
        colors = {
            (True, 0): self.SETTINGS["bands_up_linecolor"],
            (True, 1): self.SETTINGS["bands_down_linecolor"],
            (False, 0): self.SETTINGS["bands_linecolor"],
        }
        fermi_energy_mapping = {
            (False, 0): self.fermi_energy.get("fermi_energy_up", None),
            (False, 1): self.fermi_energy.get("fermi_energy_down", None),
        }

        bands_data = self.bands_data
        # Convert paths to a list of Scatter objects
        scatter_objects = []

        spin_polarized = 1 in bands_data["band_type_idx"]
        for spin in [0, 1]:
            # In case of non-spin-polarized or SOC calculations, the spin index is only 0
            if spin not in bands_data["band_type_idx"]:
                continue

            x_bands = np.array(bands_data["x"]).reshape(1, -1)
            # New shape: (number of bands, number of kpoints)
            y_bands = bands_data["y"][:, bands_data["band_type_idx"] == spin].T
            # Concatenate the bands and prepare the traces
            x_bands_comb, y_bands_comb = _prepare_combined_plotly_traces(
                x_bands, y_bands
            )

            fermi_energy = fermi_energy_mapping.get(
                ("fermi_energy" in self.fermi_energy, spin),
                self.fermi_energy.get("fermi_energy"),
            )

            scatter_objects.append(
                go.Scatter(
                    x=x_bands_comb,
                    y=y_bands_comb - fermi_energy,
                    mode="lines",
                    line={
                        "color": colors[(spin_polarized, spin)],
                        "shape": "spline",
                        "smoothing": 1.3,
                    },
                    showlegend=False,
                )
            )

        self._add_traces_to_fig(fig, scatter_objects, 1)

    def _add_pdos_traces(self, fig):
        # Extract DOS data
        dos_data = self.pdos_data["dos"]

        # Pre-allocate memory for Scatter objects
        num_traces = len(dos_data)
        scatter_objects = [None] * num_traces

        # dictionary with keys (bool(spin polarized), bool(spin up))
        fermi_energy_spin_mapping = {
            (False, True): self.fermi_energy.get("fermi_energy_up", None),
            (False, False): self.fermi_energy.get("fermi_energy_down", None),
        }

        # Vectorize Scatter object creation
        for i, trace in enumerate(dos_data):
            dos_np = np.array(trace["x"])
            fill = "tozerox" if self.plot_type == "combined" else "tozeroy"
            fermi_energy = fermi_energy_spin_mapping.get(
                ("fermi_energy" in self.fermi_energy, trace["label"].endswith("(↑)")),
                self.fermi_energy.get("fermi_energy"),
            )

            x_data = (
                trace["y"] if self.plot_type == "combined" else dos_np - fermi_energy
            )
            y_data = (
                dos_np - fermi_energy if self.plot_type == "combined" else trace["y"]
            )
            scatter_objects[i] = go.Scatter(
                x=x_data,
                y=y_data,
                fill=fill,
                name=trace["label"],
                line={
                    "color": trace["borderColor"],
                    "shape": "spline",
                    "smoothing": 1.0,
                },
                legendgroup=trace["label"],
            )

        self._add_traces_to_fig(fig, scatter_objects, 2)

    def _add_projection_traces(self, fig):
        """Function to add the projected bands traces to the bands plot."""
        projected_bands = self.bands_data["projected_bands"]
        # dictionary with keys (bool(spin polarized), bool(spin up))
        fermi_energy_spin_mapping = {
            (False, True): self.fermi_energy.get("fermi_energy_up", None),
            (False, False): self.fermi_energy.get("fermi_energy_down", None),
        }

        scatter_objects = []
        for proj_bands in projected_bands:
            fermi_energy = fermi_energy_spin_mapping.get(
                (
                    "fermi_energy" in self.fermi_energy,
                    proj_bands["label"].endswith("(↑)"),
                ),
                self.fermi_energy.get("fermi_energy"),
            )
            scatter_objects.append(
                go.Scatter(
                    x=proj_bands["x"],
                    y=np.array(proj_bands["y"]) - fermi_energy,
                    fill="toself",
                    legendgroup=proj_bands["label"],
                    mode="lines",
                    line={"width": 0, "color": proj_bands["color"]},
                    name=proj_bands["label"],
                    # If PDOS is present, use those legend entries
                    showlegend=True if self.plot_type == "bands" else False,
                )
            )

        self._add_traces_to_fig(fig, scatter_objects, 1)

    def _customize_combined_layout(self, fig):
        self._customize_layout(fig, self._bands_xaxis, self._bands_yaxis)
        self._customize_layout(fig, self._pdos_xaxis, self._pdos_yaxis, col=2)
        fig.update_layout(
            legend={"xanchor": "left", "x": 1.06},
            height=self.SETTINGS["combined_plot_height"],
            width=self.SETTINGS["combined_plot_width"],
            plot_bgcolor="white",
        )

    def _customize_layout(self, fig, xaxis, yaxis, row=1, col=1):
        fig.update_xaxes(patch=xaxis, row=row, col=col)
        fig.update_yaxes(patch=yaxis, row=row, col=col, showticklabels=True)
        fig.add_hline(
            y=0,
            line={
                "color": self.SETTINGS["horizontal_linecolor"],
                "width": 1,
                "dash": "dot",
            },
            row=row,
            col=col,
        )

    def _customize_single_layout(self, fig):
        xaxis = getattr(self, f"_{self.plot_type}_xaxis")
        yaxis = getattr(self, f"_{self.plot_type}_yaxis")

        fig.update_layout(
            xaxis=xaxis,
            yaxis=yaxis,
            plot_bgcolor="white",
            height=self.SETTINGS[f"{self.plot_type}_plot_height"],
            width=self.SETTINGS[f"{self.plot_type}_plot_width"],
        )


class BandPdosWidget(ipw.VBox):
    """
    A widget for plotting band structure and projected density of states (PDOS) data.

    Parameters:
    - bands (optional): A node containing band structure data.
    - pdos (optional): A node containing PDOS data.

    Attributes:
    - description: HTML description of the widget.
    - dos_atoms_group: Dropdown widget to select the grouping of atoms for PDOS plotting.
    - dos_plot_group: Dropdown widget to select the type of PDOS contributions to plot.
    - selected_atoms: Text widget to select specific atoms for PDOS plotting.
    - update_plot_button: Button widget to update the plot.
    - download_button: Button widget to download the data.
    - project_bands_box: Checkbox widget to choose whether projected bands should be plotted.
    - dos_data: PDOS data.
    - bands_data: Band structure data.
    - bandsplot_widget: Plotly widget for band structure and PDOS plot.
    - bands_widget: Output widget to display the bandsplot widget.
    - pdos_options_out: Output widget to clear specific widgets.
    """

    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Select the style of plotting the projected density of states.
        </div>"""
    )
    projected_bands_width = 0.5

    def __init__(self, bands=None, pdos=None, **kwargs):
        if bands is None and pdos is None:
            raise ValueError("Either bands or pdos must be provided")

        self.bands = bands  # bands node
        self.pdos = pdos  # pdos node

        self.dos_atoms_group = ipw.Dropdown(
            description="Group by:",
            options=[
                ("Kinds", "kinds"),
                ("Atoms", "atoms"),
            ],
            value="kinds",
            style={"description_width": "initial"},
        )
        self.dos_plot_group = ipw.Dropdown(
            description="Plot contributions:",
            options=[
                ("Total", "total"),
                ("Orbital", "orbital"),
                ("Angular momentum", "angular_momentum"),
            ],
            value="total",
            style={"description_width": "initial"},
        )
        self.selected_atoms = ipw.Text(
            description="Select atoms:",
            value="",
            style={"description_width": "initial"},
        )
        self._wrong_syntax = StatusHTML(clear_after=8)
        self.update_plot_button = ipw.Button(
            description="Update Plot",
            icon="pencil",
            button_style="primary",
            disabled=False,
        )
        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(visibility="hidden"),
        )
        self.project_bands_box = ipw.Checkbox(
            value=False,
            description="Add `fat bands` projections",
        )

        # Information for the plot
        self.pdos_data = self._get_pdos_data()
        self.bands_data = self._get_bands_data()
        # Plotly widget
        self.bandsplot_widget = BandPdosPlotly(
            bands_data=self.bands_data, pdos_data=self.pdos_data
        ).bandspdosfigure
        # Output widget to display the bandsplot widget
        self.bands_widget = ipw.Output()
        # Output widget to clear the specific widgets
        self.pdos_options_out = ipw.Output()

        pdos_options_list = [
            self.description,
            self.dos_atoms_group,
            self.dos_plot_group,
            ipw.HBox([self.selected_atoms, self._wrong_syntax]),
            self.update_plot_button,
        ]
        # If projections are available in the bands data, include the box to plot fat-bands
        if self.bands_data and "projected_bands" in self.bands_data:
            pdos_options_list.insert(4, self.project_bands_box)

        self.pdos_options = ipw.VBox(pdos_options_list)

        self._initial_view()

        # Set the event handlers
        self.download_button.on_click(self.download_data)
        self.update_plot_button.on_click(self._update_plot)

        super().__init__(
            children=[
                self.pdos_options_out,
                self.download_button,
                self.bands_widget,  # Add the output widget to the VBox
            ],
            **kwargs,
        )

        # Plot the options only if the pdos is provided or in case the bands data contains projections
        if self.pdos or (self.bands_data and "projected_bands" in self.bands_data):
            with self.pdos_options_out:
                display(self.pdos_options)

    def download_data(self, _=None):
        """Function to download the data."""
        file_name_bands = "bands_data.json"
        file_name_pdos = "dos_data.json"
        if self.bands_data:
            bands_data_export = {}
            for key, value in self.bands_data.items():
                if isinstance(value, np.ndarray):
                    bands_data_export[key] = value.tolist()
                else:
                    bands_data_export[key] = value

            json_str = json.dumps(bands_data_export)
            b64_str = base64.b64encode(json_str.encode()).decode()
            self._download(payload=b64_str, filename=file_name_bands)
        if self.pdos_data:
            json_str = json.dumps(self.pdos_data)
            b64_str = base64.b64encode(json_str.encode()).decode()
            self._download(payload=b64_str, filename=file_name_pdos)

    @staticmethod
    def _download(payload, filename):
        """Download payload as a file named as filename."""
        from IPython.display import Javascript

        javas = Javascript(
            f"""
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """
        )
        display(javas)

    def _get_pdos_data(self):
        if not self.pdos:
            return None
        expanded_selection, syntax_ok = string_range_to_list(
            self.selected_atoms.value, shift=-1
        )
        if syntax_ok:
            pdos = get_pdos_data(
                self.pdos,
                group_tag=self.dos_atoms_group.value,
                plot_tag=self.dos_plot_group.value,
                selected_atoms=expanded_selection,
            )
            return pdos
        return None

    def _get_bands_data(self):
        if not self.bands:
            return None

        expanded_selection, syntax_ok = string_range_to_list(
            self.selected_atoms.value, shift=-1
        )
        if syntax_ok:
            bands = get_bands_projections_data(
                self.bands,
                group_tag=self.dos_atoms_group.value,
                plot_tag=self.dos_plot_group.value,
                selected_atoms=expanded_selection,
                bands_width=self.projected_bands_width,
            )
            return bands
        return None

    def _initial_view(self):
        with self.bands_widget:
            self._clear_output_and_display(self.bandsplot_widget)
            self.download_button.layout.visibility = "visible"
            self.project_bands_box.layout.visibility = "visible"

    def _update_plot(self, _=None):
        with self.bands_widget:
            expanded_selection, syntax_ok = string_range_to_list(
                self.selected_atoms.value, shift=-1
            )
            if not syntax_ok:
                self._wrong_syntax.message = """<div class='alert alert-danger'> ERROR: Invalid syntax for selected atoms</div>"""
                clear_output(wait=True)
            else:
                self.pdos_data = self._get_pdos_data()
                self.bands_data = self._get_bands_data()
                self.bandsplot_widget = BandPdosPlotly(
                    bands_data=self.bands_data,
                    pdos_data=self.pdos_data,
                    project_bands=self.project_bands_box.value,
                ).bandspdosfigure
                self._clear_output_and_display(self.bandsplot_widget)

    def _clear_output_and_display(self, widget=None):
        clear_output(wait=True)
        if widget:
            display(widget)


def _prepare_combined_plotly_traces(x_to_conc, y_to_conc):
    """Combine multiple lines into a single trace.

    The rows of y are concatenated with a np.nan column as a separator. Moreover,
    the x values are ajduced to match the shape of the concatenated y values. These
    transfomred arrays, representing multiple datasets/lines, can be plotted in a single trace.
    """
    if y_to_conc.ndim != 2:
        raise ValueError("y must be a 2D array")

    y_dim0 = y_to_conc.shape[0]

    # Add a np.nan column as a separator
    y_transf = np.hstack(
        [
            y_to_conc,
            np.full((y_dim0, 1), np.nan),
        ]
    ).flatten()

    # Same logic for the x axis
    x_transf = x_to_conc.reshape(1, -1) * np.ones(y_dim0).reshape(-1, 1)
    x_transf = np.hstack([x_transf, np.full((y_dim0, 1), np.nan)]).flatten()

    return x_transf, y_transf


def _prepare_projections_to_plot(bands_data, projections, bands_width):
    """Prepare the projected bands to be plotted.

    This function transforms the projected bands into a format that can be plotted
    in a single trace. To use the fill option `toself`,
    a band needs to be concatenated with its mirror image, first.
    """
    projected_bands = []
    for spin in [0, 1]:
        # In case of non-spin-polarized calculations, the spin index is only 0
        if spin not in bands_data["band_type_idx"]:
            continue

        x_bands = bands_data["x"]
        # New shape: (number of bands, number of kpoints)
        y_bands = bands_data["y"][:, bands_data["band_type_idx"] == spin].T

        for proj in projections[spin]:
            # Create the upper and lower boundary of the fat bands based on the orbital projections
            y_bands_proj_upper = y_bands + bands_width * proj["projections"].T
            y_bands_proj_lower = y_bands - bands_width * proj["projections"].T
            # As mentioned above, the bands need to be concatenated with their mirror image
            # to create the filled areas properly
            y_bands_mirror = np.hstack(
                [y_bands_proj_upper, y_bands_proj_lower[:, ::-1]]
            )
            # Same logic for the energy axis
            x_bands_mirror = np.concatenate([x_bands, x_bands[::-1]]).reshape(1, -1)
            x_bands_comb, y_bands_proj_comb = _prepare_combined_plotly_traces(
                x_bands_mirror, y_bands_mirror
            )

            projected_bands.append(
                {
                    "x": x_bands_comb.tolist(),
                    "y": y_bands_proj_comb.tolist(),
                    "label": proj["label"],
                    "color": proj["color"],
                }
            )
    return projected_bands


def get_bands_projections_data(
    outputs, group_tag, plot_tag, selected_atoms, bands_width, fermi_energy=None
):
    """Extract the bandstructure and possibly the projections along the bands."""
    if "band_structure" not in outputs:
        return None

    bands_data = outputs.band_structure._get_bandplot_data(
        cartesian=True, prettify_format=None, join_symbol=None, get_segments=True
    )
    # The fermi energy from band calculation is not robust.
    if "fermi_energy_up" in outputs.band_parameters:
        bands_data["fermi_energy_up"] = outputs.band_parameters["fermi_energy_up"]
        bands_data["fermi_energy_down"] = outputs.band_parameters["fermi_energy_down"]
    else:
        bands_data["fermi_energy"] = (
            outputs.band_parameters["fermi_energy"] or fermi_energy
        )

    bands_data["pathlabels"] = get_bands_labeling(bands_data)

    if "projwfc" in outputs:
        projections = []

        if "projections" in outputs.projwfc:
            projections.append(
                _projections_curated_options(
                    outputs.projwfc.projections,
                    spin_type="none",
                    group_tag=group_tag,
                    plot_tag=plot_tag,
                    selected_atoms=selected_atoms,
                    projections_pdos="projections",
                )
            )
        else:
            for spin_proj, spin_type in zip(
                [
                    outputs.projwfc.projections_up,
                    outputs.projwfc.projections_down,
                ],
                ["up", "down"],
            ):
                projections.append(
                    _projections_curated_options(
                        spin_proj,
                        spin_type=spin_type,
                        group_tag=group_tag,
                        plot_tag=plot_tag,
                        selected_atoms=selected_atoms,
                        projections_pdos="projections",
                    )
                )

        bands_data["projected_bands"] = _prepare_projections_to_plot(
            bands_data, projections, bands_width
        )
    return bands_data


def get_pdos_data(pdos, group_tag, plot_tag, selected_atoms):
    dos = []

    if "output_dos" not in pdos.dos:
        return None

    _, energy_dos, _ = pdos.dos.output_dos.get_x()
    tdos_values = {f"{n}": v for n, v, _ in pdos.dos.output_dos.get_y()}

    if "projections" in pdos.projwfc:
        # Total DOS
        tdos = {
            "label": "Total DOS",
            "x": energy_dos.tolist(),
            "y": tdos_values.get("dos").tolist(),
            "borderColor": "#8A8A8A",  # dark gray
            "backgroundColor": "#999999",  # light gray
            "backgroundAlpha": "40%",
            "lineStyle": "solid",
        }
        dos.append(tdos)
        dos += _projections_curated_options(
            pdos.projwfc.projections,
            spin_type="none",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )
    else:
        # Total DOS (↑) and Total DOS (↓)
        tdos_up = {
            "label": "Total DOS (↑)",
            "x": energy_dos.tolist(),
            "y": tdos_values.get("dos_spin_up").tolist(),
            "borderColor": "#8A8A8A",  # dark gray
            "backgroundColor": "#999999",  # light gray
            "backgroundAlpha": "40%",
            "lineStyle": "solid",
        }
        tdos_down = {
            "label": "Total DOS (↓)",
            "x": energy_dos.tolist(),
            "y": (-tdos_values.get("dos_spin_down")).tolist(),
            "borderColor": "#8A8A8A",  # dark gray
            "backgroundColor": "#999999",  # light gray
            "backgroundAlpha": "40%",
            "lineStyle": "dash",
        }
        dos += [tdos_up, tdos_down]

        # Spin-up (↑) and Spin-down (↓)
        dos += _projections_curated_options(
            pdos.projwfc.projections_up,
            spin_type="up",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )
        dos += _projections_curated_options(
            pdos.projwfc.projections_down,
            spin_type="down",
            line_style="dash",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )

    data_dict = {
        "dos": dos,
    }
    if "fermi_energy_up" in pdos.nscf.output_parameters:
        data_dict["fermi_energy_up"] = pdos.nscf.output_parameters["fermi_energy_up"]
        data_dict["fermi_energy_down"] = pdos.nscf.output_parameters[
            "fermi_energy_down"
        ]
    else:
        data_dict["fermi_energy"] = pdos.nscf.output_parameters["fermi_energy"]

    return json.loads(json.dumps(data_dict))


def _get_grouping_key(
    group_tag,
    plot_tag,
    atom_position,
    kind_name,
    orbital_name_plotly,
    orbital_angular_momentum,
):
    """Generates the grouping key based on group_tag and plot_tag."""

    key_formats = {
        ("atoms", "total"): r"{var1}-{var}",
        ("kinds", "total"): r"{var1}",
        ("atoms", "orbital"): r"{var1}-{var}<br>{var2}",
        ("kinds", "orbital"): r"{var1}-{var2}",
        ("atoms", "angular_momentum"): r"{var1}-{var}<br>{var3}",
        ("kinds", "angular_momentum"): r"{var1}-{var3}",
    }

    key = key_formats.get((group_tag, plot_tag))
    if key is not None:
        return key.format(
            var=atom_position,
            var1=kind_name,
            var2=orbital_name_plotly,
            var3=orbital_angular_momentum,
        )
    else:
        return None


def _curate_orbitals(orbital):
    """Curate and transform the orbital data into the desired format."""
    # Constants for HTML tags
    HTML_TAGS = {
        "s": "s",
        "pz": "p<sub>z</sub>",
        "px": "p<sub>x</sub>",
        "py": "p<sub>y</sub>",
        "dz2": "d<sub>z<sup>2</sup></sub>",
        "dxy": "d<sub>xy</sub>",
        "dxz": "d<sub>xz</sub>",
        "dyz": "d<sub>yz</sub>",
        "dx2-y2": "d<sub>x<sup>2</sup>-y<sup>2</sup></sub>",
        "fz3": "f<sub>z<sup>3</sup></sub>",
        "fxz2": "f<sub>xz<sup>2</sup></sub>",
        "fyz2": "f<sub>yz<sup>2</sup></sub>",
        "fxyz": "f<sub>xzy</sub>",
        "fx(x2-3y2)": "f<sub>x(x<sup>2</sup>-3y<sup>2</sup>)</sub>",
        "fy(3x2-y2)": "f<sub>y(3x<sup>2</sup>-y<sup>2</sup>)</sub>",
        "fy(x2-z2)": "f<sub>y(x<sup>2</sup>-z<sup>2</sup>)</sub>",
        0.5: "<sup>+1</sup>/<sub>2</sub>",
        -0.5: "<sup>-1</sup>/<sub>2</sub>",
        1.5: "<sup>+3</sup>/<sub>2</sub>",
        -1.5: "<sup>-3</sup>/<sub>2</sub>",
        2.5: "<sup>+5</sup>/<sub>2</sub>",
        -2.5: "<sup>-5</sup>/<sub>2</sub>",
    }

    orbital_data = orbital.get_orbital_dict()
    kind_name = orbital_data["kind_name"]
    atom_position = [round(i, 2) for i in orbital_data["position"]]

    try:
        orbital_name = orbital.get_name_from_quantum_numbers(
            orbital_data["angular_momentum"], orbital_data["magnetic_number"]
        ).lower()
        orbital_name_plotly = HTML_TAGS.get(orbital_name, orbital_name)
        orbital_angular_momentum = orbital_name[0]
    except AttributeError:
        # Set quanutum numbers
        qn_j = orbital_data["total_angular_momentum"]
        qn_l = orbital_data["angular_momentum"]
        qn_m_j = orbital_data["magnetic_number"]
        orbital_name = f"j {qn_j} l {qn_l} m_j{qn_m_j}"
        orbital_name_plotly = f"j={HTML_TAGS.get(qn_j, qn_j)} <i>l</i>={qn_l} m<sub>j</sub>={HTML_TAGS.get(qn_m_j, qn_m_j)}"
        orbital_angular_momentum = f"l {qn_l} "

    return orbital_name_plotly, orbital_angular_momentum, kind_name, atom_position


def _projections_curated_options(
    projections: ProjectionData,
    group_tag,
    plot_tag,
    selected_atoms,
    projections_pdos="pdos",
    spin_type="none",
    line_style="solid",
):
    """Extract and curate the projections.

    This function can be used to extract the PDOS or the projections data.
    """
    _proj_pdos = {}
    list_positions = []

    # Constants for spin types
    SPIN_LABELS = {"up": "(↑)", "down": "(↓)", "none": ""}
    SIGN_MULT_FACTOR = {"up": 1, "down": -1, "none": 1}

    if projections_pdos == "pdos":
        proj_data = projections.get_pdos()
    elif projections_pdos == "projections":
        proj_data = projections.get_projections()
    else:
        raise ValueError(f"Invalid value for `projections_pdos`: {projections_pdos}")

    for orb_proj in proj_data:
        if projections_pdos == "pdos":
            orbital, proj_pdos, energy = orb_proj
        elif projections_pdos == "projections":
            orbital, proj_pdos = orb_proj
            energy = None

        (
            orbital_name_plotly,
            orbital_angular_momentum,
            kind_name,
            atom_position,
        ) = _curate_orbitals(orbital)

        if atom_position not in list_positions:
            list_positions.append(atom_position)

        key = _get_grouping_key(
            group_tag,
            plot_tag,
            atom_position,
            kind_name,
            orbital_name_plotly,
            orbital_angular_momentum,
        )
        if not selected_atoms:
            if key:
                _proj_pdos.setdefault(key, [energy, 0])[1] += proj_pdos

        else:
            try:
                index = list_positions.index(atom_position)
                if index in selected_atoms:
                    if key:
                        _proj_pdos.setdefault(key, [energy, 0])[1] += proj_pdos

            except ValueError:
                pass

    curated_proj = []
    for label, (energy, proj_pdos) in _proj_pdos.items():
        label += SPIN_LABELS[spin_type]  # noqa: PLW2901
        if projections_pdos == "pdos":
            orbital_proj_pdos = {
                "label": label,
                "x": energy.tolist(),
                "y": (SIGN_MULT_FACTOR[spin_type] * proj_pdos).tolist(),
                "borderColor": cmap(label),
                "lineStyle": line_style,
            }
        else:
            orbital_proj_pdos = {
                "label": label,
                "projections": proj_pdos,
                "color": cmap(label),
            }
        curated_proj.append(orbital_proj_pdos)

    return curated_proj


def get_bands_labeling(bandsdata: dict) -> list:
    """Function to return two lists containing the labels and values (kpoint) for plotting.
    params:
    - bandsdata: dictionary from `get_bands_projections_data` function
    output:  update bandsdata with a new key "pathlabels" including (list of str), label_values (list of float)
    """
    UNICODE_SYMBOL = {
        "GAMMA": "\u0393",
        "DELTA": "\u0394",
        "LAMBDA": "\u039b",
        "SIGMA": "\u03a3",
        "EPSILON": "\u0395",
    }
    paths = bandsdata.get("paths")
    labels = []
    for path in paths:  # Remove duplicates
        label_a = [path["from"], path["x"][0]]
        label_b = [path["to"], path["x"][-1]]
        if label_a not in labels:
            labels.append(label_a)
        if label_b not in labels:
            labels.append(label_b)

    clean_labels = []  # Format
    for i in labels:
        if clean_labels:
            if (i not in clean_labels) and (clean_labels[-1][-1] == i[1]):
                clean_labels[-1][0] = clean_labels[-1][0] + "|" + i[0]
            else:
                clean_labels.append(i)
        else:
            clean_labels.append(i)

    path_labels = [label[0] for label in clean_labels]
    for i, label in enumerate(path_labels):
        path_labels[i] = re.sub(
            r"([A-Z]+)", lambda x: UNICODE_SYMBOL.get(x.group(), x.group()), label
        )
    path_values = [label[1] for label in clean_labels]
    return [path_labels, path_values]


def cmap(label: str) -> str:
    """Return RGB string of color for given pseudo info
    Hardcoded at the momment.
    """
    import random

    # if a unknow type generate random color based on ascii sum
    ascn = sum([ord(c) for c in label])
    random.seed(ascn)

    return f"#{random.randint(0, 0xFFFFFF):06x}"
