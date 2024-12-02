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
    points, equivalent_sites_data, gamma=0.3, sigma=0.3, _label="", intensity=1.0
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
            intensity = equivalent_sites_data[site]["multiplicity"] * intensity
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
