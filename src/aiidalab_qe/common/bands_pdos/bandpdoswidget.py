import ipywidgets as ipw
from IPython.display import clear_output, display

from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base.utils import StatusHTML, string_range_to_list

from .bandpdosplotly import BandPdosPlotly
from .model import BandsPdosModel


class BandPdosWidget(ipw.VBox):
    """
    A widget for plotting band structure and projected density of states (PDOS) data.

    Parameters
    ----------
    - bands (optional): A node containing band structure data.
    - pdos (optional): A node containing PDOS data.

    Attributes
    ----------
    - description: HTML description of the widget.
    - dos_atoms_group: Dropdown widget to select the grouping of atoms for PDOS plotting.
    - dos_plot_group: Dropdown widget to select the type of PDOS contributions to plot.
    - selected_atoms: Text widget to select specific atoms for PDOS plotting.
    - update_plot_button: Button widget to update the plot.
    - download_button: Button widget to download the data.
    - project_bands_box: Checkbox widget to choose whether projected bands should be plotted.
    - plot_widget: Plotly widget for band structure and PDOS plot.
    - bands_widget: Output widget to display the bandsplot widget.
    """

    def __init__(self, model: BandsPdosModel, bands=None, pdos=None, **kwargs):
        if bands is None and pdos is None:
            raise ValueError("Either bands or pdos must be provided")

        super().__init__(
            children=[LoadingWidget("Loading widgets")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_needs_bands_projections_change,
            "needs_projections_controls",
        )
        self._model.observe(
            self._on_needs_pdos_options_change,
            "needs_pdos_options",
        )

        self.rendered = False

        self._model.bands = bands
        self._model.pdos = pdos

    def render(self):
        if self.rendered:
            return

        self.dos_atoms_group = ipw.Dropdown(
            description="Group by:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._model, "dos_atoms_group_options"),
            (self.dos_atoms_group, "options"),
        )
        ipw.link(
            (self._model, "dos_atoms_group"),
            (self.dos_atoms_group, "value"),
        )
        self.dos_atoms_group.observe(
            self._update_plot,
            "value",
        )

        self.dos_plot_group = ipw.Dropdown(
            description="Plot contributions:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._model, "dos_plot_group_options"),
            (self.dos_plot_group, "options"),
        )
        ipw.link(
            (self._model, "dos_plot_group"),
            (self.dos_plot_group, "value"),
        )
        self.dos_plot_group.observe(
            self._update_plot,
            "value",
        )

        self.selected_atoms = ipw.Text(
            placeholder="e.g. 1..5 8 10",
            description="Select atoms:",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "selected_atoms"),
            (self.selected_atoms, "value"),
        )

        self._wrong_syntax = StatusHTML(clear_after=8)

        self.update_plot_button = ipw.Button(
            description="Apply selection",
            icon="pencil",
            button_style="primary",
        )
        self.update_plot_button.on_click(self._update_plot)

        self.download_button = ipw.Button(
            description="Download Data",
            icon="download",
            button_style="primary",
            layout=ipw.Layout(visibility="hidden"),
        )
        self.download_button.on_click(self._model.download_data)

        self.project_bands_box = ipw.Checkbox(
            description="Add `fat bands` projections",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "project_bands_box"),
            (self.project_bands_box, "value"),
        )
        self.project_bands_box.observe(
            self._update_plot,
            "value",
        )

        self.proj_bands_width_slider = ipw.FloatSlider(
            min=0.01,
            max=2.0,
            step=0.01,
            description="`Fat bands` max width (eV):",
            orientation="horizontal",
            readout=True,
            readout_format=".2f",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="380px", visibility="hidden"),
        )
        ipw.link(
            (self._model, "proj_bands_width"),
            (self.proj_bands_width_slider, "value"),
        )
        self.proj_bands_width_slider.observe(
            self._update_plot,
            "value",
        )

        self.proj_controls = ipw.VBox(
            children=[
                self.project_bands_box,
                self.proj_bands_width_slider,
            ],
            layout=ipw.Layout(display="none"),
        )

        self.pdos_options = ipw.VBox(
            children=[
                ipw.HTML("""
                    <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px; max-width: 600px;">
                        Select the style of plotting the projected density of states.
                    </div>
                """),
                self.dos_atoms_group,
                self.dos_plot_group,
                ipw.HBox(
                    children=[
                        self.selected_atoms,
                        self._wrong_syntax,
                        self.update_plot_button,
                    ]
                ),
                self.proj_controls,
            ],
            layout=ipw.Layout(display="none"),
        )

        self.legend_interaction_description = ipw.HTML(
            """
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 5px; max-width: 600px;">
                    The legend entries can be clicked to hide or show the corresponding
                    data. Double-clicking on a legend entry will isolate it.
                </div>
            """,
            layout=ipw.Layout(display="none"),
        )

        # Output widget to display the bandsplot widget
        self.bands_widget = ipw.Output()

        self.children = [
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px; max-width: 600px;">
                    Hover over the plot to reveal controls for zoom, pan, and
                    downloading the image. Use the zoom tools or your mouse to zoom in
                    on specific regions, and click on the axes for interactive features.
                    The home button resets to the default view, and the autoscale option
                    displays all computed data, including semicore states.
                </div>
            """),
            self.pdos_options,
            self.download_button,
            self.legend_interaction_description,
            self.bands_widget,
        ]

        self.rendered = True

        self._initial_view()
        self._toggle_projection_controls()
        self._toggle_pdos_options()
        self._update_plot()

    def _on_needs_bands_projections_change(self, _):
        self._toggle_projection_controls()

    def _on_needs_pdos_options_change(self, _):
        self._toggle_pdos_options()

    def _initial_view(self):
        with self.bands_widget:
            self._clear_output_and_display()
            self.download_button.layout.visibility = "visible"
            self.project_bands_box.layout.visibility = "visible"

    def _toggle_projection_controls(self):
        """If projections are available in the bands data,
        show the box to plot fat-bands."""
        if not self.rendered:
            return
        self.proj_controls.layout.display = (
            "flex" if self._model.needs_projections_controls else "none"
        )

    def _toggle_pdos_options(self):
        """Plot the options only if the pdos is provided or in case the bands data
        contains projections."""
        if not self.rendered:
            return
        if self._model.needs_pdos_options:
            self.pdos_options.layout.display = "flex"
            self.legend_interaction_description.layout.display = "flex"
        else:
            self.pdos_options.layout.display = "none"
            self.legend_interaction_description.layout.display = "none"

    def _update_plot(self, _=None):
        with self.bands_widget:
            clear_output(wait=True)
            display(LoadingWidget("Plotting results"))
            _, syntax_ok = string_range_to_list(self._model.selected_atoms, shift=-1)
            if not syntax_ok:
                self._wrong_syntax.message = """
                    <div class='alert alert-danger'>
                        ERROR: Invalid syntax for selected atoms
                    </div>
                """
                clear_output(wait=True)
            else:
                self._model.fetch_data()

                self.plot = BandPdosPlotly(
                    bands_data=self._model.bands_data,
                    pdos_data=self._model.pdos_data,
                    project_bands=self._model.project_bands_box,
                ).bandspdosfigure

                # Get current axis range
                xaxis_range = list(self.plot.layout["xaxis"]["range"])
                yaxis_range = list(self.plot.layout["yaxis"]["range"])

                self._clear_output_and_display()

                # Restore Old axis range. I do it after the plot is displayed to the Reset button always return to the Default SETTINGs
                if self._model.bands_data:
                    self.plot.plotly_relayout({"yaxis.range": yaxis_range})
                if self._model.pdos_data and not self._model.bands_data:
                    self.plot.plotly_relayout({"xaxis.range": xaxis_range})

                self.proj_bands_width_slider.layout.visibility = (
                    "visible" if self._model.project_bands_box else "hidden"
                )

    def _clear_output_and_display(self):
        clear_output(wait=True)
        if hasattr(self, "plot"):
            display(self.plot)
