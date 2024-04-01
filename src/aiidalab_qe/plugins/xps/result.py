"""XPS results view widgets"""

import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultPanel


def export_xps_data(outputs):
    """Export the data from the XPS workchain"""

    chemical_shifts = {}
    symmetry_analysis_data = outputs.symmetry_analysis_data.get_dict()
    equivalent_sites_data = symmetry_analysis_data["equivalent_sites_data"]
    if "chemical_shifts" in outputs:
        for key, data in outputs.chemical_shifts.items():
            ele = key[:-4]
            chemical_shifts[ele] = data.get_dict()
    binding_energies = {}
    if "binding_energies" in outputs:
        for key, data in outputs.binding_energies.items():
            ele = key[:-3]
            binding_energies[ele] = data.get_dict()

    return (
        chemical_shifts,
        binding_energies,
        equivalent_sites_data,
    )


def xps_spectra_broadening(
    points, equivalent_sites_data, gamma=0.3, sigma=0.3, label=""
):
    """Broadening the XPS spectra with Voigt function and return the spectra data"""

    import numpy as np
    from scipy.special import voigt_profile  # pylint: disable=no-name-in-module

    result_spectra = {}
    fwhm_voight = gamma / 2 + np.sqrt(gamma**2 / 4 + sigma**2)
    for element, point in points.items():
        result_spectra[element] = {}
        final_spectra_y_arrays = []
        total_multiplicity = sum(
            [equivalent_sites_data[site]["multiplicity"] for site in point]
        )
        max_core_level_shift = max(point.values())
        min_core_level_shift = min(point.values())
        # Energy range for the Broadening function
        x_energy_range = np.linspace(
            min_core_level_shift - fwhm_voight - 1.5,
            max_core_level_shift + fwhm_voight + 1.5,
            500,
        )
        for site in point:
            # Weight for the spectra of every atom
            intensity = equivalent_sites_data[site]["multiplicity"]
            relative_core_level_position = point[site]
            y = (
                intensity
                * voigt_profile(
                    x_energy_range - relative_core_level_position, sigma, gamma
                )
                / total_multiplicity
            )
            result_spectra[element][site] = [x_energy_range, y]
            final_spectra_y_arrays.append(y)
        total = sum(final_spectra_y_arrays)
        result_spectra[element]["total"] = [x_energy_range, total]
    return result_spectra


class Result(ResultPanel):
    title = "XPS"
    workchain_labels = ["xps"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)
        self.experimental_data = None  # Placeholder for experimental data

    def _update_view(self):
        import plotly.graph_objects as go

        spectrum_select_prompt = ipw.HTML(
            """
            <div style="line-height: 140%; padding-top: 10px; padding-right: 10px; padding-bottom: 0px;"><b>
            Select spectrum to plot</b></div>"""
        )

        voigt_profile_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                Set the <a href="https://en.wikipedia.org/wiki/Voigt_profile" target="_blank">Voigt profile</a> to broaden the XPS spectra:
            </div>"""
        )

        spectra_type = ipw.ToggleButtons(
            options=[
                ("Chemical shift", "chemical_shift"),
                ("Binding energy", "binding_energy"),
            ],
            value="chemical_shift",
        )
        gamma = ipw.FloatSlider(
            value=0.1,
            min=0.01,
            max=0.5,
            description="Lorentzian profile ($\gamma$)",
            disabled=False,
            style={"description_width": "initial"},
        )
        sigma = ipw.FloatSlider(
            value=0.1,
            min=0.01,
            max=0.5,
            description="Gaussian profile ($\sigma$)",
            disabled=False,
            style={"description_width": "initial"},
        )
        fill = ipw.Checkbox(
            description="Fill",
            value=True,
            disabled=False,
            style={"description_width": "initial"},
        )
        # Create a description label
        upload_description = ipw.HTML(
            value="<b>Upload Experimental Data (csv format):</b>",
            placeholder="",
            description="",
        )

        # Create the upload button
        upload_btn = ipw.FileUpload(
            description="Choose File",
            multiple=False,
        )
        upload_container = ipw.VBox([upload_description, upload_btn])
        upload_btn.observe(self._handle_upload, names="value")

        paras = ipw.HBox(
            children=[
                gamma,
                sigma,
            ]
        )
        # get data
        (
            chemical_shifts,
            binding_energies,
            equivalent_sites_data,
        ) = export_xps_data(self.outputs.xps)
        self.spectrum_select_options = [
            key.split("_")[0] for key in chemical_shifts.keys()
        ]
        spectrum_select = ipw.Dropdown(
            description="",
            disabled=False,
            value=self.spectrum_select_options[0],
            options=self.spectrum_select_options,
            layout=ipw.Layout(width="20%"),
        )
        # init figure
        self.g = go.FigureWidget(
            layout=go.Layout(
                title=dict(text="XPS"),
                barmode="overlay",
            )
        )
        self.g.layout.xaxis.title = "Chemical shift (eV)"
        self.g.layout.xaxis.autorange = "reversed"
        #
        spectra = xps_spectra_broadening(
            chemical_shifts, equivalent_sites_data, gamma=gamma.value, sigma=sigma.value
        )
        # only plot the selected spectrum
        for site, d in spectra[spectrum_select.value].items():
            self.g.add_scatter(
                x=d[0], y=d[1], fill="tozeroy", name=site.replace("_", " ")
            )

        def response(change):
            data = []
            if spectra_type.value == "chemical_shift":
                points = chemical_shifts
                xaxis = "Chemical Shift (eV)"
            else:
                points = binding_energies
                xaxis = "Binding Energy (eV)"
            #
            spectra = xps_spectra_broadening(
                points, equivalent_sites_data, gamma=gamma.value, sigma=sigma.value
            )

            for site, d in spectra[spectrum_select.value].items():
                data.append(
                    {
                        "x": d[0],
                        "y": d[1],
                        "site": site,
                    }
                )
            fill_type = "tozeroy" if fill.value else None
            with self.g.batch_update():
                if len(self.g.data) == len(data):
                    for i in range(len(data)):
                        self.g.data[i].x = data[i]["x"]
                        self.g.data[i].y = data[i]["y"]
                        self.g.data[i].fill = fill_type
                        self.g.data[i].name = data[i]["site"].replace("_", " ")

                else:
                    self.g.data = []
                    for d in data:
                        self.g.add_scatter(
                            x=d["x"], y=d["y"], fill=fill_type, name=d["site"]
                        )
                self.g.layout.barmode = "overlay"
                self.g.layout.xaxis.title = xaxis
            self.plot_experimental_data()

        spectra_type.observe(response, names="value")
        spectrum_select.observe(response, names="value")
        gamma.observe(response, names="value")
        sigma.observe(response, names="value")
        fill.observe(response, names="value")
        self.children = [
            spectra_type,
            ipw.HBox(
                children=[
                    spectrum_select_prompt,
                    spectrum_select,
                ]
            ),
            voigt_profile_help,
            paras,
            fill,
            self.g,
            upload_container,
        ]

    def _handle_upload(self, change):
        """Process the uploaded experimental data file."""
        import pandas as pd

        uploaded_file = next(iter(change.new.values()))
        content = uploaded_file["content"]
        content_str = content.decode("utf-8")

        from io import StringIO

        df = pd.read_csv(StringIO(content_str), header=None)

        self.experimental_data = df
        self.plot_experimental_data()

    def plot_experimental_data(self):
        """Plot the experimental data alongside the calculated data."""
        if self.experimental_data is not None:
            x = self.experimental_data[0]
            y = self.experimental_data[1]
            self.g.add_scatter(x=x, y=y, mode="lines", name="Experimental Data")
