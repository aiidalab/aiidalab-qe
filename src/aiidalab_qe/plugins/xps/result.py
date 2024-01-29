"""XPS results view widgets

"""
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
    spectra_cls = {}
    if "xps_spectra_cls" in outputs:
        for key, data in outputs.xps_spectra_cls.items():
            ele = key[:-12]
            X = data.get_x()[1]
            Y = data.get_y()
            array_y = []
            for y in Y:
                array_y.append(y[1])

            # The total dos parsed
            spectra_cls[ele] = {"x": X, "y": array_y}
    spectra_be = {}
    if "xps_spectra_be" in outputs:
        for key, data in outputs.xps_spectra_be.items():
            ele = key[:-11]
            X = data.get_x()[1]
            Y = data.get_y()
            array_y = []
            for y in Y:
                array_y.append(y[1])

            # The total dos parsed
            spectra_be[ele] = {"x": X, "y": array_y}

    return (
        chemical_shifts,
        binding_energies,
        spectra_cls,
        spectra_be,
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
            relative_peak_position = point[site]
            y = (
                intensity
                * voigt_profile(x_energy_range - relative_peak_position, sigma, gamma)
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

    def _update_view(self):
        import plotly.graph_objects as go

        voigt_profile_help = ipw.HTML(
            """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        The Voigt function:
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
            value=0.3,
            min=0.1,
            max=1,
            description="γ",
            disabled=False,
            style={"description_width": "initial"},
        )
        sigma = ipw.FloatSlider(
            value=0.3,
            min=0.1,
            max=1,
            description="σ",
            disabled=False,
            style={"description_width": "initial"},
        )
        fill = ipw.Checkbox(
            description="Fill",
            value=True,
            disabled=False,
            style={"description_width": "initial"},
        )
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
            spectra_cls,
            spectra_be,
            equivalent_sites_data,
        ) = export_xps_data(self.outputs.xps)
        # init figure
        g = go.FigureWidget(
            layout=go.Layout(
                title=dict(text="XPS"),
                barmode="overlay",
            )
        )
        g.layout.xaxis.title = "Chemical shift (eV)"
        g.layout.xaxis.autorange = "reversed"
        #
        spectra = xps_spectra_broadening(
            chemical_shifts, equivalent_sites_data, gamma=gamma.value, sigma=sigma.value
        )
        for element, data in spectra.items():
            # print(element, data)
            for site, d in data.items():
                g.add_scatter(x=d[0], y=d[1], fill="tozeroy", name=f"{element}_{site}")

        def response(change):
            X = []
            Y = []
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
            for _key, data in spectra.items():
                for _site, d in data.items():
                    X.append(d[0])
                    Y.append(d[1])

            with g.batch_update():
                for i in range(len(X)):
                    g.data[i].x = X[i]
                    g.data[i].y = Y[i]
                    if fill.value:
                        g.data[i].fill = "tozeroy"
                    else:
                        g.data[i].fill = None
                g.layout.barmode = "overlay"
                g.layout.xaxis.title = xaxis

        spectra_type.observe(response, names="value")
        gamma.observe(response, names="value")
        sigma.observe(response, names="value")
        fill.observe(response, names="value")
        self.children = [
            spectra_type,
            voigt_profile_help,
            paras,
            fill,
            g,
        ]
