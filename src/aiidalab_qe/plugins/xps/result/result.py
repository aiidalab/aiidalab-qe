"""XPS results view widgets"""

import ipywidgets as ipw
import plotly.graph_objects as go

from aiidalab_qe.common.panel import ResultsPanel

from .model import XpsResultsModel


class XpsResultsPanel(ResultsPanel[XpsResultsModel]):
    experimental_data = None  # Placeholder for experimental data

    def _on_file_upload(self, change):
        self._model.upload_experimental_data(change["new"])

    def _render(self):
        spectra_type = ipw.ToggleButtons()
        ipw.dlink(
            (self._model, "spectra_type_options"),
            (spectra_type, "options"),
        )
        ipw.link(
            (self._model, "spectra_type"),
            (spectra_type, "value"),
        )
        spectra_type.observe(
            self._update_plot,
            "value",
        )

        gamma = ipw.FloatSlider(
            min=0.01,
            max=0.5,
            step=0.01,
            description=r"Lorentzian profile ($\gamma$)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "gamma"),
            (gamma, "value"),
        )
        gamma.observe(
            self._update_plot,
            "value",
        )

        sigma = ipw.FloatSlider(
            min=0.01,
            max=0.5,
            step=0.01,
            description=r"Gaussian profile ($\sigma$)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "sigma"),
            (sigma, "value"),
        )
        sigma.observe(
            self._update_plot,
            "value",
        )

        self.intensity = ipw.FloatText(
            min=0.001,
            description="Adjustable Intensity Factor",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "intensity"),
            (self.intensity, "value"),
        )
        self.intensity.observe(
            self._update_plot,
            "value",
        )

        fill = ipw.Checkbox(
            description="Fill",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "fill"),
            (fill, "value"),
        )
        fill.observe(
            self._update_plot,
            "value",
        )

        # Create a description label
        upload_description = ipw.HTML(
            value="Upload Experimental Data (<b>csv format, without header</b>):",
            placeholder="",
            description="",
        )

        # Create the upload button
        upload_btn = ipw.FileUpload(
            description="Choose File",
            multiple=False,
        )
        upload_btn.observe(
            self._on_file_upload,
            "value",
        )

        upload_container = ipw.VBox(
            children=[
                upload_description,
                upload_btn,
                self.intensity,
            ],
        )

        parameters_container = ipw.HBox(
            children=[
                gamma,
                sigma,
                fill,
            ]
        )

        self.spectrum_select = ipw.Dropdown(
            description="",
            disabled=False,
            layout=ipw.Layout(width="20%"),
        )
        ipw.dlink(
            (self._model, "spectrum_options"),
            (self.spectrum_select, "options"),
        )
        ipw.link(
            (self._model, "spectrum"),
            (self.spectrum_select, "value"),
        )
        self.spectrum_select.observe(
            self._update_plot,
            "value",
        )

        self.plot = go.FigureWidget(
            layout=go.Layout(
                title={"text": "XPS"},
                barmode="overlay",
            )
        )
        self.plot.layout.xaxis.title = "Chemical shift (eV)"
        self.plot.layout.xaxis.autorange = "reversed"

        self.results_container.children = [
            spectra_type,
            ipw.HBox(
                children=[
                    ipw.HTML("""
                        <div style="line-height: 140%; padding-top: 10px; padding-right: 10px; padding-bottom: 0px;">
                            <b>Select spectrum to plot</b>
                        </div>
                    """),
                    self.spectrum_select,
                ]
            ),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                    Set the <a href="https://en.wikipedia.org/wiki/Voigt_profile" target="_blank">Voigt profile</a> to broaden the XPS spectra:
                </div>
            """),
            parameters_container,
            self.plot,
            upload_container,
        ]

    def _post_render(self):
        self._model.update_spectrum_options()

    def _update_plot(self, _):
        if not self.rendered:
            return

        data, x_axis_label, fill_type = self._model.get_data()

        with self.plot.batch_update():
            if len(self.plot.data) == len(data):
                for i in range(len(data)):
                    self.plot.data[i].x = data[i]["x"]
                    self.plot.data[i].y = data[i]["y"]
                    self.plot.data[i].fill = fill_type
                    self.plot.data[i].name = data[i]["site"].replace("_", " ")

            else:
                self.plot.data = []
                for d in data:
                    self.plot.add_scatter(
                        x=d["x"],
                        y=d["y"],
                        fill=fill_type,
                        name=d["site"],
                    )

            self.plot.layout.barmode = "overlay"
            self.plot.layout.xaxis.title = x_axis_label

        self._plot_experimental_data()

    def _plot_experimental_data(self):
        """Plot the experimental data alongside the calculated data."""
        if not self.rendered:
            return
        if self._model.experimental_data:
            x = self._model.experimental_data[0]
            y = self._model.experimental_data[1]
            self.plot.add_scatter(x=x, y=y, mode="lines", name="Experimental Data")
