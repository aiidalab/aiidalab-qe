import ipywidgets as ipw

from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base.utils import StatusHTML, string_range_to_list

from .model import BandsPdosModel


class BandsPdosWidget(ipw.VBox):
    """A widget for plotting band structure and projected density of states (PDOS).

    Parameters
    ----------
    `model`: `BandsPdosModel`
        The MVC model containing the data and logic for the widget.

    Attributes
    ----------
    `description`: `ipywidgets.HTML`
         HTML description of the widget.
    `dos_atoms_group`: `ipywidgets.Dropdown`
         Dropdown widget to select the grouping of atoms for PDOS plotting.
    `dos_plot_group`: `ipywidgets.Dropdown`
         Dropdown widget to select the type of PDOS contributions to plot.
    `selected_atoms`: `ipywidgets.Text`
         Text widget to select specific atoms for PDOS plotting.
    `update_plot_button`: `ipywidgets.Button`
         Button widget to update the plot.
    `download_button`: `ipywidgets.Button`
         Button widget to download the data.
    `project_bands_box`: `ipywidgets.Checkbox`
         Checkbox widget to choose whether projected bands should be plotted.
    `plot`: `plotly.graph_objects.FigureWidget`
         Plotly widget for band structure and PDOS plot.
    """

    def __init__(self, model: BandsPdosModel, **kwargs):
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

    def render(self):
        if self.rendered:
            return

        self.dos_atoms_group = ipw.Dropdown(
            description="Atom grouping:",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="350px"),
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
            self._update_pdos_plot,
            "value",
        )

        self.dos_plot_group = ipw.Dropdown(
            description="Orbital grouping:",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="350px"),
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
            self._update_pdos_plot,
            "value",
        )

        self.selected_atoms = ipw.Text(
            placeholder="e.g. 1..5 8 10",
            description="Select atoms:",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="350px"),
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
        self.update_plot_button.on_click(self._update_pdos_plot)

        self.download_button = ipw.Button(
            description="Download data",
            icon="download",
            button_style="primary",
            layout=ipw.Layout(visibility="hidden"),
        )
        self.download_button.on_click(self._model.download_data)

        self.download_image = ipw.Button(
            description="Download image",
            button_style="primary",
            icon="fa-image",
        )
        self.download_image.on_click(self._model.download_image)
        self.image_format = ipw.Dropdown(
            description="Format:",
            layout=ipw.Layout(width="auto"),
        )
        ipw.dlink((self._model, "image_format_options"), (self.image_format, "options"))
        ipw.link((self._model, "image_format"), (self.image_format, "value"))

        self.download_buttons = ipw.HBox(
            children=[self.download_button, self.download_image, self.image_format]
        )
        self.project_bands_box = ipw.Checkbox(
            description="Add `fat bands` projections",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "project_bands_box"),
            (self.project_bands_box, "value"),
        )
        self.project_bands_box.observe(
            self._update_bands_projections,
            "value",
        )

        self.proj_bands_width_slider = ipw.FloatSlider(
            min=0.01,
            max=5.0,
            step=0.01,
            description="`Fat bands` max width (eV):",
            orientation="horizontal",
            continuous_update=False,
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
            self._update_bands_projections_thickness,
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

        self.color_picker = ipw.ColorPicker(description="Pick a color")
        self.trace_selector = ipw.Dropdown(description="Select a trace")

        self.color_selector = ipw.VBox(
            children=[
                ipw.HTML("""
                    <div style="line-height: 140%; padding: 10px 0; max-width: 600px;">
                        <strong>How to customize the plot:</strong><br>
                        Select a trace from the dropdown menu, adjust its color using the color picker, and see the plot update instantly.
                    </div>
                    """),
                ipw.HBox(
                    children=[
                        self.trace_selector,
                        self.color_picker,
                    ]
                ),
            ],
        )
        ipw.dlink(
            (self._model, "trace_selector_options"),
            (self.trace_selector, "options"),
        )
        ipw.link(
            (self._model, "color_picker"),
            (self.color_picker, "value"),
        )
        self.trace_selector.observe(
            self._trace_selector_change,
            "value",
        )
        self.color_picker.observe(
            self._update_trace_color,
            "value",
        )

        # Aspect ratio
        self.horizontal_width_percentage = ipw.IntSlider(
            min=30,
            max=100,
            step=5,
            description="Horizonal width %:",
            orientation="horizontal",
            continuous_update=False,
            readout=True,
            readout_format=".0f",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="380px"),
        )
        ipw.link(
            (self._model, "horizontal_width_percentage"),
            (self.horizontal_width_percentage, "value"),
        )
        self.horizontal_width_percentage.observe(
            self._on_horizontal_width_change,
            "value",
        )

        self.bands_width_percentage = ipw.IntSlider(
            min=10,
            max=90,
            step=5,
            description="Bands width %:",
            orientation="horizontal",
            continuous_update=False,
            readout=True,
            readout_format=".0f",
            style={"description_width": "initial"},
            layout=ipw.Layout(width="380px"),
        )
        ipw.link(
            (self._model, "bands_width_percentage"),
            (self.bands_width_percentage, "value"),
        )
        self.bands_width_percentage.observe(
            self._on_bands_width_change,
            "value",
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
            self.download_buttons,
            self.legend_interaction_description,
        ]

        self.rendered = True

        self._initial_plot()
        self._toggle_projection_controls()
        self._toggle_pdos_options()

    def _on_needs_bands_projections_change(self, _):
        self._toggle_projection_controls()

    def _on_needs_pdos_options_change(self, _):
        self._toggle_pdos_options()

    def _initial_plot(self):
        """Create the initial plot."""
        self._model.fetch_data()
        self._model.create_plot()
        self.plot = self._model.plot
        self.proj_bands_width_slider.layout.visibility = (
            "visible" if self._model.project_bands_box else "hidden"
        )
        self.download_button.layout.visibility = "visible"
        self.project_bands_box.layout.visibility = "visible"
        self.children = [
            *self.children,
            ipw.Box(
                children=[self.plot],
                layout=ipw.Layout(margin="0 auto"),
            ),
            self.color_selector,
            self.horizontal_width_percentage,
        ]
        if self._model.helper.plot_type == "combined":
            self.children = [
                *self.children,
                self.bands_width_percentage,
            ]

    def _update_bands_projections(self, _):
        """Update the plot with the selected projection."""
        self.proj_bands_width_slider.layout.visibility = (
            "visible" if self._model.project_bands_box else "hidden"
        )

        _, syntax_ok = string_range_to_list(self._model.selected_atoms, shift=-1)
        if not syntax_ok:
            self._wrong_syntax.message = """
                <div class='alert alert-danger'>
                    ERROR: Invalid syntax for selected atoms
                </div>
            """
        else:
            self._model.update_bands_projections()
            self._trace_selector_change({"new": 0})

    def _update_bands_projections_thickness(self, _):
        """Update the plot with the selected projection thickness."""
        self._model.update_bands_projections_thickness()

    def _update_pdos_plot(self, _):
        """Update the plot with the selected PDOS options."""
        _, syntax_ok = string_range_to_list(self._model.selected_atoms, shift=-1)
        if not syntax_ok:
            self._wrong_syntax.message = """
                <div class='alert alert-danger'>
                    ERROR: Invalid syntax for selected atoms
                </div>
            """
        else:
            self._model.update_pdos_plot()
            self._trace_selector_change({"new": 0})

    def _toggle_projection_controls(self):
        """If projections are available in the bands data,
        show the box to plot fat-bands."""
        if not self.rendered:
            return
        self.proj_controls.layout.display = (
            "flex" if self._model.needs_projections_controls else "none"
        )
        self._model.project_bands_box = True

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

    def _trace_selector_change(self, change):
        self._model.update_color_picker(change["new"])

    def _update_trace_color(self, change):
        self._model.update_trace_color(change["new"])

    def _on_horizontal_width_change(self, change):
        self._model.update_horizontal_width(change["new"])

    def _on_bands_width_change(self, change):
        self._model.update_column_width_ratio(change["new"])
