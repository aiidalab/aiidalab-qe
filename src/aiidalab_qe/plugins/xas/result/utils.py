import numpy as np
from scipy.interpolate import make_interp_spline


def write_csv(dataset):
    from pandas import DataFrame

    x_vals = dataset[0]["x"]
    df_data = {"energy_ev": x_vals}
    for entry in dataset:
        if "site" in entry["name"]:
            if entry["weighting"] != 1:
                df_data[
                    f'{entry["name"].capitalize().replace("_", " ")} (Weighted)'
                ] = entry["y"]
                df_data[
                    f'{entry["name"].capitalize().replace("_", " ")} (Unweighted)'
                ] = entry["y"] / entry["weighting"]
            else:
                df_data[entry["name"].capitalize().replace("_", " ")] = entry["y"]
        else:
            df_data[entry["name"]] = entry["y"]

    df = DataFrame(data=df_data)
    df_energy_indexed = df.set_index("energy_ev")

    return df_energy_indexed.to_csv(header=True)


def export_xas_data(outputs):
    final_spectra = {}
    equivalent_sites_data = {}
    if "final_spectra" in outputs:
        final_spectra = outputs.final_spectra
        symmetry_analysis_data = outputs.symmetry_analysis_data.get_dict()
        equivalent_sites_data = symmetry_analysis_data["equivalent_sites_data"]
    return (
        final_spectra,
        equivalent_sites_data,
    )


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

            if y <= 1.0e-6:  # do this to skip the calculation for very small values
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
