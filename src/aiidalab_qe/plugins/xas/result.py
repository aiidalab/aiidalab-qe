"""XAS results view widgets

"""
import base64
import hashlib
from typing import Callable

import ipywidgets as ipw
import numpy as np
from IPython.display import HTML, display
from scipy.interpolate import make_interp_spline

from aiidalab_qe.common.panel import ResultPanel


class SpectrumDownloadButton(ipw.Button):
    """Download button with dynamic content
    The content is generated using a callback when the button is clicked.
    Modified from responses to https://stackoverflow.com/questions/61708701/how-to-download-a-file-using-ipywidget-button#62641240
    """

    def __init__(self, filename: str, contents: Callable[[], str], **kwargs):
        super(SpectrumDownloadButton, self).__init__(**kwargs)
        self.filename = filename
        self.contents = contents
        self.on_click(self.__on_click)

    def __on_click(self, b):
        if self.contents() is None:
            pass  # to avoid a crash because NoneType obviously can't be processed here
        else:
            contents: bytes = self.contents().encode("utf-8")
            b64 = base64.b64encode(contents)
            payload = b64.decode()
            digest = hashlib.md5(contents).hexdigest()  # bypass browser cache
            id = f"dl_{digest}"

            display(
                HTML(
                    f"""
                <html>
                <body>
                <a id="{id}" download="{self.filename}" href="data:text/csv;base64,{payload}" download>
                </a>
                <script>
                (function download() {{
                document.getElementById('{id}').click();
                }})()
                </script>
                </body>
                </html>
            """
                )
            )


def write_csv(dataset):
    from pandas import DataFrame

    x_vals = dataset[0]["x"]
    df_data = {"energy_ev": x_vals}
    # header = []
    for entry in dataset:
        # df_data[f'{entry["name"]}(weight_{entry["weighting_string"]})'] = entry["y"]
        if "site" in entry["name"]:
            if entry["weighting"] != 1:
                df_data[
                    f'{entry["name"].capitalize().replace("_", " ")} (Weighted)'
                ] = entry["y"]
                df_data[
                    f'{entry["name"].capitalize().replace("_", " ")} (Unweighted)'
                ] = (entry["y"] / entry["weighting"])
            else:
                df_data[entry["name"].capitalize().replace("_", " ")] = entry["y"]
        else:
            df_data[entry["name"]] = entry["y"]
        # header.append(f'weighting={str(entry["weighting_string"])}')

    df = DataFrame(data=df_data)
    df_energy_indexed = df.set_index("energy_ev")

    return df_energy_indexed.to_csv(header=True)


def export_xas_data(outputs):
    if "final_spectra" in outputs.xas:
        final_spectra = outputs.xas.final_spectra
        symmetry_analysis_data = outputs.xas.symmetry_analysis_data.get_dict()
        equivalent_sites_data = symmetry_analysis_data["equivalent_sites_data"]

        return (
            final_spectra,
            equivalent_sites_data,
        )
    else:
        return None


def broaden_xas(
    input_array, variable=False, gamma_hole=0.01, gamma_max=5, center_energy=15
):
    """Take an input spectrum and return a broadened spectrum as output using either a constant or variable parameter.

    :param input_array: The 2D array of x/y values to be broadened. Should be plotted with
                        little or no broadening before using the function.
    :param gamma_hole: The broadening parameter for the Lorenzian broadening function. In constant mode (variable=False),
                       this value is applied to the entire spectrum. In variable mode (variable=True), this value defines
                       the starting broadening parameter of the arctangent function. Refers to the natural linewidth of
                       the element/XAS edge combination in question and (for elements Z > 10) should be based on reference
                       values from X-ray spectroscopy.
    :param variable: Request a variable-energy broadening of the spectrum instead of the defaultconstant broadening.
                     Uses the functional form defined in Calandra and Bunau, PRB, 87, 205105 (2013).
    :param gamma_max: The maximum lorenzian broadening to be applied in variable energy broadening mode. Refers to the
                      broadening effects at infinite energy above the main edge.
    :param center_energy: The inflection point of the variable broadening function. Does not relate to experimental data
                          and must be tuned manually.
    """

    if variable:
        if not all([gamma_hole, gamma_max, center_energy]):
            missing = [
                i[0]
                for i in zip(
                    ["gamma_hole", "gamma_max", "center_energy"],
                    [gamma_hole, gamma_max, center_energy],
                )
                if i[1] is None
            ]
            raise ValueError(
                f"The following variables were not defined {missing} and are required for variable-energy broadening"
            )

    x_vals = input_array[:, 0]
    y_vals = input_array[:, 1]
    # y_vals_norm = np.array(y_vals/np.trapz(y_vals, x_vals))

    lorenz_y = np.zeros(len(x_vals))

    if variable:
        for x, y in zip(x_vals, y_vals):
            if x < 0:  # the function is bounded between gamma_hole and gamma_max
                gamma_var = gamma_hole
            else:
                e = x / center_energy

                gamma_var = gamma_hole + gamma_max * (
                    0.5 + np.arctan((e - 1) / (e**2)) / np.pi
                )

            if x <= 1.0e-6:  # do this to avoid trying to broaden values close to 0
                lorenz_y = y
            else:
                lorenz_y += (
                    gamma_var
                    / 2.0
                    / np.pi
                    / ((x_vals - x) ** 2 + 0.25 * gamma_var**2)
                    * y
                )
    else:
        for x, y in zip(x_vals, y_vals):
            lorenz_y += (
                gamma_hole
                / 2.0
                / np.pi
                / ((x_vals - x) ** 2 + 0.25 * gamma_hole**2)
                * y
            )

    return np.column_stack((x_vals, lorenz_y))


def get_aligned_spectra(core_wc_dict, equivalent_sites_dict):
    """Return a set of spectra aligned according to the chemical shift (difference in Fermi level).

    Primarily this is a copy of ``get_spectra_by_element`` from AiiDA-QE which operates on only one
    element.
    """
    data_dict = {}
    spectrum_dict = {
        site: node.outputs.powder_spectrum for site, node in core_wc_dict.items()
    }
    for key, value in core_wc_dict.items():
        xspectra_out_params = value.outputs.parameters_xspectra__xas_0.get_dict()
        energy_zero = xspectra_out_params["energy_zero"]
        multiplicity = equivalent_sites_dict[key]["multiplicity"]

        if "total_multiplicity" not in data_dict:
            data_dict["total_multiplicity"] = multiplicity
        else:
            data_dict["total_multiplicity"] += multiplicity

        data_dict[key] = {
            "spectrum_node": spectrum_dict[key],
            "multiplicity": multiplicity,
            "energy_zero": energy_zero,
        }

    spectra_list = []
    total_multiplicity = data_dict.pop("total_multiplicity")
    for key in data_dict:
        spectrum_node = data_dict[key]["spectrum_node"]
        site_multiplicity = data_dict[key]["multiplicity"]
        weighting = site_multiplicity / total_multiplicity
        weighting_string = f"{site_multiplicity}/{total_multiplicity}"
        spectrum_x = spectrum_node.get_x()[1]
        spectrum_y = spectrum_node.get_y()[0][1]
        spline = make_interp_spline(spectrum_x, spectrum_y)
        norm_y = spline(spectrum_x) / np.trapz(spline(spectrum_x), spectrum_x)
        weighted_spectrum = np.column_stack(
            (spectrum_x, norm_y * (site_multiplicity / total_multiplicity))
        )
        spectra_list.append(
            (
                weighted_spectrum,
                key,
                weighting,
                weighting_string,
                float(data_dict[key]["energy_zero"]),
            )
        )

    # Sort according to Fermi level, then correct to align all spectra to the
    # highest value. Note that this is needed because XSpectra automatically aligns the
    # final spectrum such that the system's Fermi level is at 0 eV.
    spectra_list.sort(key=lambda entry: entry[-1])
    highest_level = spectra_list[0][-1]
    energy_zero_corrections = [
        (entry[0], entry[1], entry[2], entry[3], entry[-1] - highest_level)
        for entry in spectra_list
    ]
    aligned_spectra = [
        (
            entry[1],
            entry[2],
            entry[3],
            np.column_stack((entry[0][:, 0] - entry[-1], entry[0][:, 1])),
        )
        for entry in energy_zero_corrections
    ]

    return aligned_spectra


class Result(ResultPanel):
    title = "XAS"
    workchain_labels = ["xas"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, identifier="xas", **kwargs)

    def _update_view(self):
        import plotly.graph_objects as go

        gamma_select_prompt = ipw.HTML(
            """
                    <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;"><b>
                    Select parameters for spectrum broadening </b></div>"""
        )

        # PNOG: If someone knows a way to format certain words differently in HTML without causing a line-break, hit me up.
        # For now, (17/10/23) we'll just have the function terms be in italics.
        # Alternatively, if it's possible to format mathematical terms in HTML without a line-break, let me know
        variable_broad_select_help = ipw.HTML(
            """
            <div style="line-height: 140%; padding-top: 5px; padding-bottom: 5px; opacity:0.5;"><b>
            Broadening parameters:
            </b>
                <p style="padding-bottom: 5px"> <i style="font-size:14pt;font-family:Times New Roman">&Gamma;<sub><em>hole</em></sub></i> - Defines a constant Lorenzian broadening width for the whole spectrum. In "variable" mode, defines the initial broadening width of the ArcTangent function. </p>
                <p style="padding-bottom: 5px"> <i style="font-size:14pt;font-family:Times New Roman">&Gamma;<sub><em>max</em></sub></i> - Maximum Lorenzian broadening parameter at infinte energy in "variable" mode.</p>
                <p style="padding-bottom: 5px"> <i style="font-size:14pt;font-family:Times New Roman"><emph>E<sub>max</em></sub></i> - Defines the inflection point of the variable-energy broadening function.</p>
            </div>
            <div style="line-height: 140%; padding-top: 5px; padding-bottom: 5px; opacity:0.5;"><b>
            Note that setting <i style="font-size:14pt;font-family:Times New Roman">&Gamma;<sub><em>hole</em></sub></i> to 0 eV will simply plot the raw spectrum.
            </div>
            """
        )
        spectrum_select_prompt = ipw.HTML(
            """
            <div style="line-height: 140%; padding-top: 0px; padding-bottom: 0px; opacity:0.5;"><b>
            Select spectrum to plot</b></div>"""
        )
        final_spectra, equivalent_sites_data = export_xas_data(self.outputs)
        xas_wc = [
            n for n in self.node.called if n.process_label == "XspectraCrystalWorkChain"
        ][0]
        core_wcs = {
            n.get_metadata_inputs()["metadata"]["call_link_label"]: n
            for n in xas_wc.called
            if n.process_label == "XspectraCoreWorkChain"
        }
        core_wc_dict = {
            key.replace("_xspectra", ""): value for key, value in core_wcs.items()
        }

        spectrum_select_options = [key.split("_")[0] for key in final_spectra.keys()]

        spectrum_select = ipw.Dropdown(
            description="",
            disabled=False,
            value=spectrum_select_options[0],
            options=spectrum_select_options,
            layout=ipw.Layout(width="20%"),
        )

        variable_broad_select = ipw.Checkbox(
            value=False,
            disabled=False,
            description="Use variable energy broadening.",
            style={"description_width": "initial", "opacity": 0.5},
        )

        gamma_hole_select = ipw.FloatSlider(
            value=0.0,
            min=0.0,
            max=5,
            step=0.1,
            description="$\Gamma_{hole}$",  # noqa: W605
            disabled=False,
            continuous_update=False,
            orientation="horizontal",
            readout=True,
        )

        gamma_max_select = ipw.FloatSlider(
            value=5.0,
            min=2.0,
            max=10,
            step=0.5,
            continuous_update=False,
            description="$\Gamma_{max}$",  # noqa: W605
            disabled=True,
            orientation="horizontal",
            readout=True,
        )

        center_e_select = ipw.FloatSlider(
            value=15.0,
            min=5,
            max=30,
            step=0.5,
            continuous_update=False,
            description="$E_{center}$",
            disabled=True,
            orientation="horizontal",
            readout=True,
        )
        download_data = SpectrumDownloadButton(
            filename="spectra.csv",
            contents=None,
            description="Download CSV",
            icon="download",
        )
        #     # get data
        #     # init figure
        g = go.FigureWidget(
            layout=go.Layout(
                title=dict(text="XAS"),
                barmode="overlay",
            )
        )

        g.layout.xaxis.title = "Relative Photon Energy (eV)"

        chosen_spectrum = spectrum_select.value
        chosen_spectrum_label = f"{chosen_spectrum}_xas"
        spectra = final_spectra[chosen_spectrum_label]

        raw_spectrum = np.column_stack((spectra.get_x()[1], spectra.get_y()[0][1]))

        x = raw_spectrum[:, 0]
        y = raw_spectrum[:, 1]
        spline = make_interp_spline(x, y)
        norm_y = spline(x) / np.trapz(spline(x), x)
        element = chosen_spectrum_label.split("_")[0]
        element_sites = [
            key
            for key in equivalent_sites_data
            if equivalent_sites_data[key]["symbol"] == element
        ]
        element_core_wcs = {}
        total_multiplicity = 0
        for site in element_sites:
            site_multiplicity = equivalent_sites_data[site]["multiplicity"]
            total_multiplicity += site_multiplicity
            element_core_wcs[site] = core_wc_dict[site]

        g.add_scatter(x=x, y=norm_y, name=f"{element} K-edge")
        for entry in get_aligned_spectra(
            core_wc_dict=element_core_wcs, equivalent_sites_dict=equivalent_sites_data
        ):
            g.add_scatter(
                x=entry[-1][:, 0],
                y=entry[-1][:, 1],
                name=entry[0].capitalize().replace("_", " "),
            )

        def _update_download_selection(dataset, element):
            download_data.contents = lambda: write_csv(dataset)
            download_data.filename = f"{element}_XAS_Spectra.csv"

        def response(change):
            chosen_spectrum = spectrum_select.value
            chosen_spectrum_label = f"{chosen_spectrum}_xas"
            element_sites = [
                key
                for key in equivalent_sites_data
                if equivalent_sites_data[key]["symbol"] == chosen_spectrum
            ]
            element_core_wcs = {
                key: value
                for key, value in core_wc_dict.items()
                if key in element_sites
            }
            spectra = []
            final_spectrum_node = final_spectra[chosen_spectrum_label]
            final_spectrum = np.column_stack(
                (final_spectrum_node.get_x()[1], final_spectrum_node.get_y()[0][1])
            )
            final_x_vals = final_spectrum[:, 0]
            final_y_vals = final_spectrum[:, 1]
            final_spectrum_spline = make_interp_spline(final_x_vals, final_y_vals)
            final_norm_y = final_spectrum_spline(final_x_vals) / np.trapz(
                final_spectrum_spline(final_x_vals), final_x_vals
            )
            spectra.append(
                (
                    f"{chosen_spectrum} K-edge",
                    1,
                    "1",
                    np.column_stack((final_x_vals, final_norm_y)),
                )
            )
            datasets = []
            for entry in get_aligned_spectra(
                core_wc_dict=element_core_wcs,
                equivalent_sites_dict=equivalent_sites_data,
            ):
                spectra.append(entry)

            for entry in spectra:
                label = entry[0]
                weighting = entry[1]
                weighting_string = entry[2]
                raw_spectrum = entry[-1]
                x = raw_spectrum[:, 0]
                y = raw_spectrum[:, 1]
                if not variable_broad_select:
                    gamma_max_select.disabled = True
                    center_e_select.disabled = True
                else:
                    gamma_max_select.disabled = False
                    center_e_select.disabled = False

                if gamma_hole_select.value == 0.0:
                    x = raw_spectrum[:, 0]
                    y = raw_spectrum[:, 1]
                else:
                    broad_spectrum = broaden_xas(
                        raw_spectrum,
                        gamma_hole=gamma_hole_select.value,
                        gamma_max=gamma_max_select.value,
                        center_energy=center_e_select.value,
                        variable=variable_broad_select.value,
                    )
                    x = broad_spectrum[:, 0]
                    y = broad_spectrum[:, 1]

                final_spline = make_interp_spline(x, y)
                final_y_vals = final_spline(final_x_vals)
                datasets.append(
                    {
                        "x": final_x_vals,
                        "y": final_y_vals,
                        "name": label,
                        "weighting": weighting,
                        "weighting_string": weighting_string,
                    }
                )
            _update_download_selection(datasets, chosen_spectrum)

            with g.batch_update():
                # If the number of datasets is different from one update to the next,
                # then we need to reset the data already in the Widget. Otherwise, we can
                # simply override the data. This also helps since then changing the
                # broadening is much smoother.
                if len(datasets) == len(
                    g.data
                ):  # if the number of entries is the same, just update
                    for index, entry in enumerate(datasets):
                        g.data[index].x = entry["x"]
                        g.data[index].y = entry["y"]
                        if "site_" in entry["name"]:
                            g.data[index].name = (
                                entry["name"].capitalize().replace("_", " ")
                            )
                        else:
                            g.data[index].name = entry["name"]
                else:  # otherwise, reset the figure
                    g.data = ()
                    for entry in datasets:
                        if "site_" in entry["name"]:
                            name = entry["name"].capitalize().replace("_", " ")
                        else:
                            name = entry["name"]
                        g.add_scatter(x=entry["x"], y=entry["y"], name=name)

        spectrum_select.observe(response, names="value")
        gamma_hole_select.observe(response, names="value")
        gamma_max_select.observe(response, names="value")
        center_e_select.observe(response, names="value")
        variable_broad_select.observe(response, names="value")
        download_data.observe(response, names=["contents", "filename"])
        self.children = [
            ipw.HBox(
                [
                    ipw.VBox(
                        [
                            spectrum_select_prompt,
                            spectrum_select,
                            gamma_select_prompt,
                            gamma_hole_select,
                            gamma_max_select,
                            center_e_select,
                        ],
                        layout=ipw.Layout(width="40%"),
                    ),
                    ipw.VBox(
                        [
                            variable_broad_select,
                            variable_broad_select_help,
                        ],
                        layout=ipw.Layout(width="60%"),
                    ),
                ]
            ),
            download_data,
            g,
        ]
