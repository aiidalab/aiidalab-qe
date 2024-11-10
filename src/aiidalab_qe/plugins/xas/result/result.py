"""XAS results view widgets"""

import ipywidgets as ipw
import plotly.graph_objects as go

from aiidalab_qe.common.panel import ResultsPanel

from .model import XasResultModel
from .spectrum_button import SpectrumDownloadButton
from .utils import write_csv


class XasResult(ResultsPanel[XasResultModel]):
    title = "XAS"
    identifier = "xas"
    workchain_labels = ["xas"]

    def render(self):
        if self.rendered:
            return

        variable_broad_select = ipw.Checkbox(
            description="Use variable energy broadening.",
            style={"description_width": "initial", "opacity": 0.5},
        )
        ipw.link(
            (self._model, "variable_broad"),
            (variable_broad_select, "value"),
        )
        variable_broad_select.observe(
            self._update_plot,
            "value",
        )

        gamma_hole_select = ipw.FloatSlider(
            min=0.0,
            max=5,
            step=0.1,
            description=r"$\Gamma_{hole}$",
            continuous_update=True,
            orientation="horizontal",
            readout=True,
        )
        ipw.link(
            (self._model, "gamma_hole"),
            (gamma_hole_select, "value"),
        )
        gamma_hole_select.observe(
            self._update_plot,
            "value",
        )

        gamma_max_select = ipw.FloatSlider(
            min=2.0,
            max=10,
            step=0.5,
            continuous_update=True,
            description=r"$\Gamma_{max}$",
            orientation="horizontal",
            readout=True,
        )
        ipw.link(
            (self._model, "gamma_max"),
            (gamma_max_select, "value"),
        )
        ipw.dlink(
            (self._model, "variable_broad"),
            (gamma_max_select, "disabled"),
            lambda value: not value,
        )
        gamma_max_select.observe(
            self._update_plot,
            "value",
        )

        center_e_select = ipw.FloatSlider(
            min=5,
            max=30,
            step=0.5,
            continuous_update=True,
            description=r"$E_{center}$",
            orientation="horizontal",
            readout=True,
        )
        ipw.link(
            (self._model, "center_e"),
            (center_e_select, "value"),
        )
        ipw.dlink(
            (self._model, "variable_broad"),
            (center_e_select, "disabled"),
            lambda value: not value,
        )
        center_e_select.observe(
            self._update_plot,
            "value",
        )

        spectrum_select = ipw.Dropdown(
            description="",
            layout=ipw.Layout(width="20%"),
        )
        ipw.dlink(
            (self._model, "spectrum_options"),
            (spectrum_select, "options"),
        )
        ipw.link(
            (self._model, "spectrum"),
            (spectrum_select, "value"),
        )
        spectrum_select.observe(
            self._update_plot,
            "value",
        )

        self.download_data = SpectrumDownloadButton(
            filename=f"{spectrum_select.value}_XAS_Spectra.csv",
            contents=None,
            description="Download CSV",
            icon="download",
        )
        self.download_data.observe(
            self._update_plot,
            ["contents", "filename"],
        )

        self.plot = go.FigureWidget(
            layout=go.Layout(
                title={"text": "XAS"},
                barmode="overlay",
            )
        )
        self.plot.layout.xaxis.title = "Relative Photon Energy (eV)"

        self.children = [
            ipw.HBox(
                [
                    ipw.VBox(
                        [
                            ipw.HTML("""
                                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 0px; opacity:0.5;">
                                    <b>Select spectrum to plot</b>
                                </div>
                            """),
                            spectrum_select,
                            ipw.HTML("""
                                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;">
                                    <b>Select parameters for spectrum broadening </b>
                                </div>
                            """),
                            gamma_hole_select,
                            gamma_max_select,
                            center_e_select,
                        ],
                        layout=ipw.Layout(width="40%"),
                    ),
                    ipw.VBox(
                        [
                            variable_broad_select,
                            # PNOG: If someone knows a way to format certain words differently in HTML without causing a line-break, hit me up.
                            # For now, (17/10/23) we'll just have the function terms be in italics.
                            # Alternatively, if it's possible to format mathematical terms in HTML without a line-break, let me know
                            ipw.HTML("""
                                <div style="line-height: 140%; padding-top: 5px; padding-bottom: 5px; opacity:0.5;">
                                    <b>Broadening parameters:</b>
                                        <p style="padding-bottom: 5px">
                                            <i style="font-size:14pt;font-family:Times New Roman">&Gamma;<sub><em>hole</em></sub></i> - Defines a constant Lorenzian broadening width for the whole spectrum. In "variable" mode, defines the initial broadening width of the ArcTangent function. </p>
                                        <p style="padding-bottom: 5px"> <i style="font-size:14pt;font-family:Times New Roman">&Gamma;<sub><em>max</em></sub></i> - Maximum Lorenzian broadening parameter at infinite energy in "variable" mode.</p>
                                        <p style="padding-bottom: 5px"> <i style="font-size:14pt;font-family:Times New Roman"><em>E<sub>center</em></sub></i> - Defines the inflection point of the variable-energy broadening function.</p>
                                    </div>
                                    <div style="line-height: 140%; padding-top: 5px; padding-bottom: 5px; opacity:0.5;">
                                    <b>Note that setting <i style="font-size:14pt;font-family:Times New Roman">&Gamma;<sub><em>hole</em></sub></i> to 0 eV will simply plot the raw spectrum.</b>
                                </div>
                            """),
                        ],
                        layout=ipw.Layout(width="60%"),
                    ),
                ]
            ),
            self.download_data,
            self.plot,
        ]

        self.rendered = True

        self._model.update_spectrum_options()

    def _update_plot(self, _):
        if not self.rendered:
            return

        datasets = self._model.get_data()
        self._update_download_selection(datasets, self._model.spectrum)

        with self.plot.batch_update():
            # If the number of datasets is different from one update to the next,
            # then we need to reset the data already in the Widget. Otherwise, we can
            # simply override the data. This also helps since then changing the
            # broadening is much smoother.

            # if the number of entries is the same, just update
            if len(datasets) == len(self.plot.data):
                for index, entry in enumerate(datasets):
                    self.plot.data[index].x = entry["x"]
                    self.plot.data[index].y = entry["y"]
                    if "site_" in entry["name"]:
                        self.plot.data[index].name = (
                            entry["name"].capitalize().replace("_", " ")
                        )
                    else:
                        self.plot.data[index].name = entry["name"]

            # otherwise, reset the figure
            else:
                self.plot.data = ()
                for entry in datasets:
                    if "site_" in entry["name"]:
                        name = entry["name"].capitalize().replace("_", " ")
                    else:
                        name = entry["name"]
                    self.plot.add_scatter(x=entry["x"], y=entry["y"], name=name)

    def _update_download_selection(self, dataset, element):
        self.download_data.contents = lambda: write_csv(dataset)
        self.download_data.filename = f"{element}_XAS_Spectra.csv"
