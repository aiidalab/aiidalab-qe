import base64
import json

import ipywidgets as ipw
import numpy as np
import plotly.graph_objects as go
from aiida.orm import ProjectionData
from aiidalab_widgets_base.utils import string_range_to_list, StatusHTML
from IPython.display import clear_output, display
from plotly.subplots import make_subplots
import re


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

    def __init__(self, bands_data=None, pdos_data=None):
        self.bands_data = bands_data
        self.pdos_data = pdos_data
        self.fermi_energy = self._get_fermi_energy()

        # Plotly Axis
        # Plotly settings
        self._bands_xaxis = self._band_xaxis()
        self._bands_yaxis = self._band_yaxis()
        self._dos_xaxis = self._dos_xaxis()
        self._dos_yaxis = self._dos_yaxis()

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
            title=dict(text="Electronic Bands (eV)", standoff=1),
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

    def _dos_xaxis(self):
        """Function to return the xaxis for the dos plot."""

        if not self.pdos_data:
            return None

        if self.bands_data:
            dosxaxis = go.layout.XAxis(
                title="Density of states",
                side="bottom",
                showgrid=True,
                showline=True,
                linecolor=self.SETTINGS["axis_linecolor"],
                mirror="ticks",
                ticks="inside",
                linewidth=2,
                tickwidth=2,
                automargin=True,
            )

        else:
            dosxaxis = go.layout.XAxis(
                title="Density of states (eV)",
                showgrid=True,
                showline=True,
                linecolor=self.SETTINGS["axis_linecolor"],
                mirror="ticks",
                ticks="inside",
                linewidth=2,
                tickwidth=2,
                range=self.SETTINGS["horizontal_range_pdos"],
            )

        return dosxaxis

    def _dos_yaxis(self):
        """Function to return the yaxis for the dos plot."""

        if not self.pdos_data:
            return None

        if self.bands_data:
            dosyaxis = go.layout.YAxis(
                # title= {"text":"Density of states (eV)", "standoff": 1},
                showgrid=True,
                showline=True,
                side="right",
                mirror="ticks",
                ticks="inside",
                linewidth=2,
                tickwidth=2,
                linecolor=self.SETTINGS["axis_linecolor"],
                zerolinewidth=2,
            )

        else:
            dosyaxis = go.layout.YAxis(
                # title="Density of states (eV)",
                showgrid=True,
                showline=True,
                side="left",
                mirror="ticks",
                ticks="inside",
                linewidth=2,
                tickwidth=2,
                linecolor=self.SETTINGS["axis_linecolor"],
                zerolinewidth=2,
            )

        return dosyaxis

    def _get_bandspdos_plot(self):
        """Function to return the bands plot widget."""
        conditions = {
            (True, False): self._create_bands_only_plot,
            (False, True): self._create_dos_only_plot,
            (True, True): self._create_combined_plot,
        }

        return conditions.get((bool(self.bands_data), bool(self.pdos_data)), None)()

    def _create_bands_only_plot(self):
        """Function to return the bands plot widget."""

        fig = go.Figure()
        paths = self.bands_data.get("paths")

        self._add_band_traces(fig, paths, "bands_only")

        band_labels = self.bands_data.get("pathlabels")
        for i in band_labels[1]:
            fig.add_vline(
                x=i, line=dict(color=self.SETTINGS["vertical_linecolor"], width=1)
            )
        fig.update_layout(
            xaxis=self._bands_xaxis,
            yaxis=self._bands_yaxis,
            plot_bgcolor="white",
            height=self.SETTINGS["bands_plot_height"],
            width=self.SETTINGS["bands_plot_width"],
        )
        return go.FigureWidget(fig)

    def _create_dos_only_plot(self):
        """Function to return the pdos plot widget."""

        fig = go.Figure()
        # Extract DOS data
        self._add_dos_traces(fig, plot_type="dos_only")
        # Add a vertical line at zero energy
        fig.add_vline(
            x=0,
            line=dict(color=self.SETTINGS["vertical_linecolor"], width=1, dash="dot"),
        )

        # Update the layout of the Figure
        fig.update_layout(
            xaxis=self._dos_xaxis,
            yaxis=self._dos_yaxis,
            plot_bgcolor="white",
            height=self.SETTINGS["pdos_plot_height"],
            width=self.SETTINGS["pdos_plot_width"],
        )

        return go.FigureWidget(fig)

    def _create_combined_plot(self):
        fig = make_subplots(
            rows=1,
            cols=2,
            shared_yaxes=True,
            column_widths=self.SETTINGS["combined_column_widths"],
            horizontal_spacing=0.015,
        )
        paths = self.bands_data.get("paths")
        self._add_band_traces(fig, paths, plot_type="combined")
        self._add_dos_traces(fig, plot_type="combined")
        band_labels = self.bands_data.get("pathlabels")
        for i in band_labels[1]:
            fig.add_vline(
                x=i,
                line=dict(color=self.SETTINGS["vertical_linecolor"], width=1),
                row=1,
                col=1,
            )
        self._customize_combined_layout(fig)
        return go.FigureWidget(fig)

    def _add_band_traces(self, fig, paths, plot_type):
        paths = self.bands_data.get("paths")

        # Spin condition: True if spin-polarized False if not
        spin_type = paths[0].get("two_band_types")
        # Convert paths to a list of Scatter objects
        scatter_objects = []

        for band in paths:
            if not spin_type:
                # Non-spin-polarized case
                for bands in band["values"]:
                    bands_np = np.array(bands)
                    scatter_objects.append(
                        go.Scatter(
                            x=band["x"],
                            y=bands_np - self.fermi_energy["fermi_energy"],
                            mode="lines",
                            line=dict(
                                color=self.SETTINGS["bands_linecolor"],
                                shape="spline",
                                smoothing=1.3,
                            ),
                            showlegend=False,
                        )
                    )
            else:
                half_len = len(band["values"]) // 2
                first_half = band["values"][:half_len]
                second_half = band["values"][half_len:]

                # Red line for the Spin up
                color_first_half = self.SETTINGS["bands_up_linecolor"]
                # Blue line for the Spin down
                color_second_half = self.SETTINGS["bands_down_linecolor"]
                if "fermi_energy" in self.fermi_energy:
                    for bands, color in zip(
                        (first_half, second_half), (color_first_half, color_second_half)
                    ):
                        bands_np = np.array(bands)
                        scatter_objects.append(
                            go.Scatter(
                                x=band["x"],
                                y=bands_np - self.fermi_energy["fermi_energy"],
                                mode="lines",
                                line=dict(
                                    color=color,
                                    shape="spline",
                                    smoothing=1.3,
                                ),
                                showlegend=False,
                            )
                        )
                else:
                    for bands, color, fermi_energy in zip(
                        (first_half, second_half),
                        (color_first_half, color_second_half),
                        (
                            self.fermi_energy["fermi_energy_up"],
                            self.fermi_energy["fermi_energy_down"],
                        ),
                    ):
                        for band_values in bands:
                            bands_np = np.array(band_values)
                            scatter_objects.append(
                                go.Scatter(
                                    x=band["x"],
                                    y=bands_np - fermi_energy,
                                    mode="lines",
                                    line=dict(
                                        color=color,
                                        shape="spline",
                                        smoothing=1.3,
                                    ),
                                    showlegend=False,
                                )
                            )

        if plot_type == "bands_only":
            fig.add_traces(scatter_objects)
        else:
            rows = [1] * len(scatter_objects)
            cols = [1] * len(scatter_objects)
            fig.add_traces(scatter_objects, rows=rows, cols=cols)

    def _add_dos_traces(self, fig, plot_type):
        # Extract DOS data
        dos_data = self.pdos_data["dos"]

        # Pre-allocate memory for Scatter objects
        num_traces = len(dos_data)
        scatter_objects = [None] * num_traces

        # Vectorize Scatter object creation
        for i, trace in enumerate(dos_data):
            dos_np = np.array(trace["x"])
            fill = "tozerox" if plot_type == "combined" else "tozeroy"

            if "fermi_energy" in self.fermi_energy:
                y_data = (
                    dos_np - self.fermi_energy["fermi_energy"]
                    if plot_type == "combined"
                    else trace["y"]
                )
                x_data = (
                    trace["y"]
                    if plot_type == "combined"
                    else dos_np - self.fermi_energy["fermi_energy"]
                )
            else:
                if trace["label"].endswith("(↑)"):
                    y_data = (
                        dos_np - self.fermi_energy["fermi_energy_up"]
                        if plot_type == "combined"
                        else trace["y"]
                    )
                    x_data = (
                        trace["y"]
                        if plot_type == "combined"
                        else dos_np - self.fermi_energy["fermi_energy_up"]
                    )
                else:
                    y_data = (
                        dos_np - self.fermi_energy["fermi_energy_down"]
                        if plot_type == "combined"
                        else trace["y"]
                    )
                    x_data = (
                        trace["y"]
                        if plot_type == "combined"
                        else dos_np - self.fermi_energy["fermi_energy_down"]
                    )
            scatter_objects[i] = go.Scatter(
                x=x_data,
                y=y_data,
                fill=fill,
                name=trace["label"],
                line=dict(color=trace["borderColor"], shape="spline", smoothing=1.0),
            )
        if plot_type == "dos_only":
            fig.add_traces(scatter_objects)
        else:
            rows = [1] * len(scatter_objects)
            cols = [2] * len(scatter_objects)
            fig.add_traces(scatter_objects, rows=rows, cols=cols)

    def _customize_combined_layout(self, fig):
        self._customize_layout(fig, self._bands_xaxis, self._bands_yaxis)
        self._customize_layout(fig, self._dos_xaxis, self._dos_yaxis, col=2)
        fig.update_layout(
            legend=dict(xanchor="left", x=1.06),
            height=self.SETTINGS["combined_plot_height"],
            width=self.SETTINGS["combined_plot_width"],
            plot_bgcolor="white",
        )

    def _customize_layout(self, fig, xaxis, yaxis, row=1, col=1):
        fig.update_xaxes(patch=xaxis, row=row, col=col)
        fig.update_yaxes(patch=yaxis, row=row, col=col, showticklabels=True)
        fig.add_hline(
            y=0,
            line=dict(color=self.SETTINGS["horizontal_linecolor"], width=1, dash="dot"),
            row=row,
            col=col,
        )

    @property
    def bandspdosfigure(self):
        return self._get_bandspdos_plot()


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

        # Information for the plot
        self.dos_data = self._get_dos_data()
        self.bands_data = self._get_bands_data()
        # Plotly widget
        self.bandsplot_widget = BandPdosPlotly(
            bands_data=self.bands_data, pdos_data=self.dos_data
        ).bandspdosfigure
        # Output widget to display the bandsplot widget
        self.bands_widget = ipw.Output()
        # Output widget to clear the specific widgets
        self.pdos_options_out = ipw.Output()

        self.pdos_options = ipw.VBox(
            [
                self.description,
                self.dos_atoms_group,
                self.dos_plot_group,
                ipw.HBox([self.selected_atoms, self._wrong_syntax]),
                self.update_plot_button,
            ]
        )

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
        if self.pdos:
            with self.pdos_options_out:
                display(self.pdos_options)

    def download_data(self, _=None):
        """Function to download the data."""
        file_name_bands = "bands_data.json"
        file_name_dos = "dos_data.json"
        if self.bands_data:
            json_str = json.dumps(self.bands_data)
            b64_str = base64.b64encode(json_str.encode()).decode()
            self._download(payload=b64_str, filename=file_name_bands)
        if self.dos_data:
            json_str = json.dumps(self.dos_data)
            b64_str = base64.b64encode(json_str.encode()).decode()
            self._download(payload=b64_str, filename=file_name_dos)

    @staticmethod
    def _download(payload, filename):
        """Download payload as a file named as filename."""
        from IPython.display import Javascript

        javas = Javascript(
            """
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """.format(payload=payload, filename=filename)
        )
        display(javas)

    def _get_dos_data(self):
        if not self.pdos:
            return None
        expanded_selection, syntax_ok = string_range_to_list(
            self.selected_atoms.value, shift=-1
        )
        if syntax_ok:
            dos = get_pdos_data(
                self.pdos,
                group_tag=self.dos_atoms_group.value,
                plot_tag=self.dos_plot_group.value,
                selected_atoms=expanded_selection,
            )
            return dos
        else:
            return None

    def _get_bands_data(self):
        if not self.bands:
            return None

        bands = export_bands_data(self.bands)
        return bands

    def _initial_view(self):
        with self.bands_widget:
            self._clear_output_and_display(self.bandsplot_widget)
            self.download_button.layout.visibility = "visible"

    def _update_plot(self, _=None):
        with self.bands_widget:
            expanded_selection, syntax_ok = string_range_to_list(
                self.selected_atoms.value, shift=-1
            )
            if not syntax_ok:
                self._wrong_syntax.message = """<div class='alert alert-danger'> ERROR: Invalid syntax for selected atoms</div>"""
                clear_output(wait=True)
            else:
                self.dos_data = self._get_dos_data()
                self.bandsplot_widget = BandPdosPlotly(
                    bands_data=self.bands_data, pdos_data=self.dos_data
                ).bandspdosfigure
                self._clear_output_and_display(self.bandsplot_widget)

    def _clear_output_and_display(self, widget=None):
        clear_output(wait=True)
        if widget:
            display(widget)


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


def _projections_curated_options(
    projections: ProjectionData,
    group_tag,
    plot_tag,
    selected_atoms,
    spin_type="none",
    line_style="solid",
):
    _pdos = {}
    list_positions = []

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

    # Constants for spin types
    SPIN_LABELS = {"up": "(↑)", "down": "(↓)", "none": ""}

    def get_key(
        group_tag,
        plot_tag,
        atom_position,
        kind_name,
        orbital_name_plotly,
        orbital_angular_momentum,
    ):
        """Generates the key based on group_tag and plot_tag."""

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

    for orbital, pdos, energy in projections.get_pdos():
        orbital_data = orbital.get_orbital_dict()
        kind_name = orbital_data["kind_name"]
        atom_position = [round(i, 2) for i in orbital_data["position"]]

        if atom_position not in list_positions:
            list_positions.append(atom_position)

        try:
            orbital_name = orbital.get_name_from_quantum_numbers(
                orbital_data["angular_momentum"], orbital_data["magnetic_number"]
            ).lower()
            orbital_name_plotly = HTML_TAGS.get(orbital_name, orbital_name)
            orbital_angular_momentum = orbital_name[0]
        except AttributeError:
            orbital_name = "j {j} l {l} m_j{m_j}".format(
                j=orbital_data["total_angular_momentum"],
                l=orbital_data["angular_momentum"],
                m_j=orbital_data["magnetic_number"],
            )
            orbital_name_plotly = "j={j} <i>l</i>={l} m<sub>j</sub>={m_j}".format(
                j=HTML_TAGS.get(
                    orbital_data["total_angular_momentum"],
                    orbital_data["total_angular_momentum"],
                ),
                l=orbital_data["angular_momentum"],
                m_j=HTML_TAGS.get(
                    orbital_data["magnetic_number"], orbital_data["magnetic_number"]
                ),
            )
            orbital_angular_momentum = "l {l} ".format(
                l=orbital_data["angular_momentum"],
            )

        if not selected_atoms:
            key = get_key(
                group_tag,
                plot_tag,
                atom_position,
                kind_name,
                orbital_name_plotly,
                orbital_angular_momentum,
            )

            if key:
                _pdos.setdefault(key, [energy, 0])[1] += pdos

        else:
            try:
                index = list_positions.index(atom_position)
                if index in selected_atoms:
                    key = get_key(
                        group_tag,
                        plot_tag,
                        atom_position,
                        kind_name,
                        orbital_name_plotly,
                        orbital_angular_momentum,
                    )

                    if key:
                        _pdos.setdefault(key, [energy, 0])[1] += pdos

            except ValueError:
                pass

    dos = []
    for label, (energy, pdos) in _pdos.items():
        if spin_type == "down":
            pdos = -pdos
            label += SPIN_LABELS[spin_type]

        if spin_type == "up":
            label += SPIN_LABELS[spin_type]

        orbital_pdos = {
            "label": label,
            "x": energy.tolist(),
            "y": pdos.tolist(),
            "borderColor": cmap(label),
            "lineStyle": line_style,
        }
        dos.append(orbital_pdos)

    return dos


def export_bands_data(outputs, fermi_energy=None):
    if "band_structure" not in outputs:
        return None

    data = json.loads(outputs.band_structure._exportcontent("json", comments=False)[0])
    # The fermi energy from band calculation is not robust.
    if "fermi_energy_up" in outputs.band_parameters:
        data["fermi_energy_up"] = outputs.band_parameters["fermi_energy_up"]
        data["fermi_energy_down"] = outputs.band_parameters["fermi_energy_down"]
    else:
        data["fermi_energy"] = outputs.band_parameters["fermi_energy"] or fermi_energy
    data["pathlabels"] = get_bands_labeling(data)
    return data


def get_bands_labeling(bandsdata: dict) -> list:
    """Function to return two lists containing the labels and values (kpoint) for plotting.
    params:
    - bandsdata: dictionary from export_bands_data function
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

    return "#%06x" % random.randint(0, 0xFFFFFF)
