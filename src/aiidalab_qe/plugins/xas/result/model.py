import numpy as np
import traitlets as tl
from scipy.interpolate import make_interp_spline

from aiidalab_qe.common.panel import ResultsModel

from .utils import broaden_xas, export_xas_data, get_aligned_spectra


class XasResultsModel(ResultsModel):
    identifier = "xas"

    spectrum_options = tl.List(
        trait=tl.Unicode(),
        default_value=[],
    )
    spectrum = tl.Unicode(None, allow_none=True)
    variable_broad = tl.Bool(False)
    gamma_hole = tl.Float(0.0)
    gamma_max = tl.Float(5.0)
    center_e = tl.Float(15.0)

    final_spectra = {}
    equivalent_sites_data = {}
    core_workchain_dict = {}

    total_multiplicity = 0

    _this_process_label = "XspectraCrystalWorkChain"

    def update_spectrum_options(self):
        if not (process_node := self.fetch_process_node()):
            return
        outputs = self._get_child_outputs()
        (
            self.final_spectra,
            self.equivalent_sites_data,
        ) = export_xas_data(outputs)
        xas_workchain = next(
            node
            for node in process_node.called
            if node.process_label == self._this_process_label
        )
        core_workchains = {
            node.get_metadata_inputs()["metadata"]["call_link_label"]: node
            for node in xas_workchain.called
            if node.process_label == "XspectraCoreWorkChain"
        }
        self.core_workchain_dict = {
            key.replace("_xspectra", ""): value
            for key, value in core_workchains.items()
        }
        options = [key.split("_")[0] for key in self.final_spectra.keys()]
        self.spectrum_options = options
        self.spectrum = options[0] if options else None

    def get_data(self):
        chosen_spectrum_label = f"{self.spectrum}_xas"

        element_sites = [
            key
            for key in self.equivalent_sites_data
            if self.equivalent_sites_data[key]["symbol"] == self.spectrum
        ]

        element_core_workchains = {}
        self.total_multiplicity = 0
        for site in element_sites:
            self.total_multiplicity += self.equivalent_sites_data[site]["multiplicity"]
            element_core_workchains[site] = self.core_workchain_dict[site]

        final_spectrum_node = self.final_spectra[chosen_spectrum_label]
        final_spectrum = np.column_stack(
            (
                final_spectrum_node.get_x()[1],
                final_spectrum_node.get_y()[0][1],
            )
        )
        final_x_vals = final_spectrum[:, 0]
        final_y_vals = final_spectrum[:, 1]
        final_spectrum_spline = make_interp_spline(final_x_vals, final_y_vals)
        integral = np.trapz(final_spectrum_spline(final_x_vals), final_x_vals)
        final_norm_y = final_spectrum_spline(final_x_vals) / integral

        spectra = [
            (
                f"{self.spectrum} K-edge",
                1,
                "1",
                np.column_stack((final_x_vals, final_norm_y)),
            )
        ]

        for entry in get_aligned_spectra(
            core_wc_dict=element_core_workchains,
            equivalent_sites_dict=self.equivalent_sites_data,
        ):
            spectra.append(entry)

        datasets = []
        for entry in spectra:
            label = entry[0]
            weighting = entry[1]
            weighting_string = entry[2]
            raw_spectrum = entry[-1]
            x = raw_spectrum[:, 0]
            y = raw_spectrum[:, 1]

            if self.gamma_hole == 0.0:
                x = raw_spectrum[:, 0]
                y = raw_spectrum[:, 1]
            else:
                broad_spectrum = broaden_xas(
                    raw_spectrum,
                    gamma_hole=self.gamma_hole,
                    gamma_max=self.gamma_max,
                    center_energy=self.center_e,
                    variable=self.variable_broad,
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

        return datasets
