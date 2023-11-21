"""XAS results view widgets

"""
import ipywidgets as ipw
import numpy as np
from scipy.interpolate import make_interp_spline

from aiidalab_qe.common.panel import ResultPanel


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
        final_spectra, _ = export_xas_data(self.outputs)

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

        # for spectrum_label, data in spectra.items():
        element = chosen_spectrum_label.split("_")[0]
        raw_spectrum = np.column_stack((spectra.get_x()[1], spectra.get_y()[0][1]))
        x = raw_spectrum[:, 0]
        y = raw_spectrum[:, 1]
        spline = make_interp_spline(x, y)
        norm_y = spline(x) / np.trapz(spline(x), x)
        g.add_scatter(x=x, y=norm_y, name=element)

        def response(change):
            chosen_spectrum = spectrum_select.value
            chosen_spectrum_label = f"{chosen_spectrum}_xas"
            spectra = final_spectra[chosen_spectrum_label]
            raw_spectrum = np.column_stack((spectra.get_x()[1], spectra.get_y()[0][1]))
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

            spline = make_interp_spline(x, y)
            norm_y = spline(x) / np.trapz(spline(x), x)

            g.update(data=[{"x": x, "y": norm_y, "name": chosen_spectrum}])

        spectrum_select.observe(response, names="value")
        gamma_hole_select.observe(response, names="value")
        gamma_max_select.observe(response, names="value")
        center_e_select.observe(response, names="value")
        variable_broad_select.observe(response, names="value")

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
                            # variable_broad_select_prompt,
                        ],
                        layout=ipw.Layout(width="60%"),
                    ),
                    # ipw.VBox([
                    # ], layout=ipw.Layout(width="60%")),
                ]
            ),
            g,
        ]
