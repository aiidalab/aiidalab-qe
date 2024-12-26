import ipywidgets as ipw

from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base.utils import StatusHTML, string_range_to_list

from .model import BandsPdosModel


class BandsPdosWidget(ipw.VBox):
    """A widget for plotting band structure and projected density of states (PDOS)."""

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
            self.download_button,
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
