import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from aiidalab_qe.common.bands_pdos.utils import (
    find_max_in_range,
    find_max_up_and_down,
    prepare_combined_plotly_traces,
)


class BandsPdosPlotly:
    SETTINGS = {
        "axis_linecolor": "#111111",
        "bands_linecolor": "#111111",
        "bands_up_linecolor": "rgba(205, 0, 0, 0.4)",  # Red Opacitiy 40%
        "bands_down_linecolor": "rgba(72, 118, 255, 0.4)",  # Blue Opacitiy 40%
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

    def __init__(
        self,
        bands_data=None,
        external_bands_data=None,
        pdos_data=None,
        bands_projections_data=None,
        plot_settings=None,
    ):
        self.bands_data = bands_data
        self.external_bands_data = external_bands_data
        self.pdos_data = pdos_data
        self.project_bands = bands_projections_data
        self.plot_settings = plot_settings or {}

        self.fermi_energy = self._get_fermi_energy()

        # Plotly axis objects
        self._bands_xaxis = self._band_xaxis()
        self._bands_yaxis = self._band_yaxis()
        self._pdos_xaxis = self._pdos_xaxis()
        self._pdos_yaxis = self._pdos_yaxis()

        # Plotly figure object
        self.figure_object = self._get_bandspdos_plot()

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
        return self.figure_object

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
        bandxaxis = go.layout.XAxis(
            title="k-points",
            range=[0, paths[-1]["x"][-1]],
            showgrid=True,
            showline=True,
            tickmode="array",
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

        self.adding_bands_traces(fig)
        self.adding_pdos_traces(fig)
        self.adding_projected_bands(fig)

        if self.plot_type == "combined":
            self._customize_combined_layout(fig)
        else:
            self._customize_single_layout(fig)

        return go.FigureWidget(fig)

    def adding_bands_traces(self, fig):
        if self.bands_data:
            bands_trace_settings = self.plot_settings.get("bands_trace_settings", None)
            self._add_band_traces(
                fig, bands_data=self.bands_data, trace_settings=bands_trace_settings
            )

            band_labels = self.bands_data.get("pathlabels")
            for label in band_labels[1]:
                fig.add_vline(
                    x=label,
                    line={"color": self.SETTINGS["vertical_linecolor"], "width": 1},
                )
            if self.external_bands_data:
                for key, bands_data in self.external_bands_data.items():
                    trace_settings = bands_data.get("trace_settings", {})
                    trace_settings.setdefault("name", key)
                    self._add_band_traces(
                        fig,
                        bands_data=bands_data,
                        trace_settings=trace_settings,
                    )

    def adding_pdos_traces(self, fig):
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

    def adding_projected_bands(self, fig):
        if self.project_bands:
            self._add_projection_traces(fig)

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

    def _add_band_traces(self, fig, bands_data, trace_settings=None):
        """Generate the band traces and add them to the figure."""
        from copy import deepcopy

        trace_settings = trace_settings or {}
        name = trace_settings.pop("name", "Bands")
        trace_settings.setdefault("dash", "solid")
        trace_settings.setdefault("shape", "linear")
        colors = {
            (True, 0): self.SETTINGS["bands_up_linecolor"],
            (True, 1): self.SETTINGS["bands_down_linecolor"],
            (False, 0): self.SETTINGS["bands_linecolor"],
        }
        fermi_energy_mapping = {
            (False, 0): self.fermi_energy.get("fermi_energy_up", None),
            (False, 1): self.fermi_energy.get("fermi_energy_down", None),
        }

        trace_name_mapping = {
            (False, 0): f"{name}",  # Base case: non-spin-polarized
            (True, 0): f"{name} (↑)",  # Spin-up case
            (True, 1): f"{name} (↓)",  # Spin-down case
        }

        # Convert paths to a list of Scatter objects
        scatter_objects = []

        spin_polarized = 1 in bands_data["band_type_idx"]
        for spin in [0, 1]:
            trace_settings_spin = deepcopy(trace_settings)
            # In case of non-spin-polarized or SOC calculations, the spin index is only 0
            if spin not in bands_data["band_type_idx"]:
                continue

            x_bands = np.array(bands_data["x"]).reshape(1, -1)
            # New shape: (number of bands, number of kpoints)
            y_bands = bands_data["y"][:, bands_data["band_type_idx"] == spin].T
            # Concatenate the bands and prepare the traces
            x_bands_comb, y_bands_comb = prepare_combined_plotly_traces(
                x_bands, y_bands
            )

            fermi_energy = fermi_energy_mapping.get(
                ("fermi_energy" in self.fermi_energy, spin),
                self.fermi_energy.get("fermi_energy"),
            )
            trace_settings_spin.setdefault("color", colors[(spin_polarized, spin)])
            scatter_objects.append(
                go.Scattergl(
                    x=x_bands_comb,
                    y=y_bands_comb - fermi_energy,
                    mode="lines",
                    line=trace_settings_spin,
                    showlegend=spin_polarized,
                    name=trace_name_mapping[(spin_polarized, spin)],
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
            scatter_objects[i] = go.Scattergl(  # type: ignore
                x=x_data,
                y=y_data,
                fill=fill,
                name=trace["label"],
                line={
                    "color": trace["borderColor"],
                    "shape": "linear",
                },
                legendgroup=trace["label"],
            )

        self._add_traces_to_fig(fig, scatter_objects, 2)

    def _prepare_bands_projection_traces_data(self, projected_bands):
        """
        Prepares data for adding or updating projected band traces.

        Args:
            projected_bands (list[dict]): List of projected bands data.

        Returns:
            list[dict]: List of dictionaries containing x, y (with Fermi offset), color, and label.
        """
        spin_mapping = {
            (False, True): self.fermi_energy.get("fermi_energy_up", None),
            (False, False): self.fermi_energy.get("fermi_energy_down", None),
        }

        prepared_data = []
        for proj_bands in projected_bands:
            fermi_energy = spin_mapping.get(
                (
                    "fermi_energy" in self.fermi_energy,
                    proj_bands["label"].endswith("(↑)"),
                ),
                self.fermi_energy.get("fermi_energy"),
            )
            prepared_data.append(
                {
                    "x": proj_bands["x"],
                    "y": np.array(proj_bands["y"]) - fermi_energy,
                    "color": proj_bands["color"],
                    "label": proj_bands["label"],
                }
            )
        return prepared_data

    def _add_projection_traces(self, fig):
        """Function to add the projected bands traces to the bands plot."""
        prepared_data = self._prepare_bands_projection_traces_data(self.project_bands)

        scatter_objects = [
            go.Scattergl(
                x=data["x"],
                y=data["y"],
                fill="toself",
                legendgroup=data["label"],
                mode="lines",
                line={"width": 0, "color": data["color"]},
                name=data["label"],
                # If PDOS is present, use those legend entries
                showlegend=True if self.plot_type == "bands" else False,
            )
            for data in prepared_data
        ]

        self._add_traces_to_fig(fig, scatter_objects, 1)

    def update_projected_bands_thickness(self, fig):
        """Update the projected bands thickness."""
        prepared_data = self._prepare_bands_projection_traces_data(self.project_bands)

        # Create a mapping from labels to y-data for efficient updates
        y_data_by_label = {data["label"]: data["y"] for data in prepared_data}

        with fig.batch_update():
            for trace in fig.data:
                if trace.xaxis == "x" and trace.legendgroup in y_data_by_label:
                    trace.y = y_data_by_label[trace.legendgroup]

    def _customize_combined_layout(self, fig):
        self._customize_layout(fig, self._bands_xaxis, self._bands_yaxis)
        self._customize_layout(fig, self._pdos_xaxis, self._pdos_yaxis, col=2)
        fig.update_layout(
            legend={"xanchor": "left", "x": 1.06},
            height=self.SETTINGS["combined_plot_height"],
            width=self.SETTINGS["combined_plot_width"],
            plot_bgcolor="white",
        )
        self._update_dos_layout(fig)

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
        self._update_dos_layout(fig)

    def _update_dos_layout(self, fig):
        def update_layout_spin_polarized(
            x_data_up,
            y_data_up,
            x_data_down,
            y_data_down,
            x_min,
            x_max,
            update_func,
            layout_type,
        ):
            most_negative_down, max_up = find_max_up_and_down(
                x_data_up, y_data_up, x_data_down, y_data_down, x_min, x_max
            )
            if layout_type == "layout":
                update_func(yaxis={"range": [most_negative_down * 1.10, max_up * 1.10]})
            elif layout_type == "xaxes":
                update_func(
                    patch={"range": [most_negative_down * 1.10, max_up * 1.10]},
                    row=1,
                    col=2,
                )

        def update_layout_non_spin_polarized(
            total_dos_xdata, total_dos_ydata, x_min, x_max, update_func, layout_type
        ):
            max_value = find_max_in_range(
                total_dos_xdata, total_dos_ydata, x_min, x_max
            )
            if layout_type == "layout":
                update_func(yaxis={"range": [0, max_value * 1.10]})
            elif layout_type == "xaxes":
                update_func(patch={"range": [0, max_value * 1.10]}, row=1, col=2)

        def get_x_min_max(fermi_energy):
            return (
                self.SETTINGS["horizontal_range_pdos"][0] + fermi_energy,
                self.SETTINGS["horizontal_range_pdos"][1] + fermi_energy,
            )

        def handle_spin_polarization(fermi_energy, update_func, layout_type):
            spin_polarized = "(↑)" in self.pdos_data["dos"][0]["label"]
            if not spin_polarized:
                total_dos_xdata = self.pdos_data["dos"][0]["x"]
                total_dos_ydata = self.pdos_data["dos"][0]["y"]
                x_min, x_max = get_x_min_max(fermi_energy)
                update_layout_non_spin_polarized(
                    total_dos_xdata,
                    total_dos_ydata,
                    x_min,
                    x_max,
                    update_func,
                    layout_type,
                )
            else:
                x_data_up = self.pdos_data["dos"][0]["x"]
                y_data_up = self.pdos_data["dos"][0]["y"]
                x_data_down = self.pdos_data["dos"][1]["x"]
                y_data_down = self.pdos_data["dos"][1]["y"]
                x_min, x_max = get_x_min_max(fermi_energy)
                update_layout_spin_polarized(
                    x_data_up,
                    y_data_up,
                    x_data_down,
                    y_data_down,
                    x_min,
                    x_max,
                    update_func,
                    layout_type,
                )

        # PDOS plot type
        if self.plot_type == "pdos":
            fermi_energy = self.fermi_energy.get(
                "fermi_energy", self.fermi_energy.get("fermi_energy_up")
            )
            handle_spin_polarization(fermi_energy, fig.update_layout, "layout")

        # Combined plot type
        if self.plot_type == "combined":
            fermi_energy = self.fermi_energy.get(
                "fermi_energy", self.fermi_energy.get("fermi_energy_up")
            )
            handle_spin_polarization(fermi_energy, fig.update_xaxes, "xaxes")
