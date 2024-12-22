import traitlets as tl

from aiidalab_qe.common.panel import ResultsModel

from .utils import export_xps_data, xps_spectra_broadening


class XpsResultsModel(ResultsModel):
    title = "XPS"
    identifier = "xps"

    spectra_type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["Chemical shift", "chemical_shift"],
            ["Binding energy", "binding_energy"],
        ],
    )
    spectra_type = tl.Unicode("chemical_shift")
    gamma = tl.Float(0.1)
    sigma = tl.Float(0.1)
    intensity = tl.Float(1.0)
    fill = tl.Bool(True)
    spectrum_options = tl.List(
        trait=tl.Unicode(),
        default_value=[],
    )
    spectrum = tl.Unicode(None, allow_none=True)

    spectra = {}
    chemical_shifts = {}
    binding_energies = {}
    equivalent_sites_data = {}
    experimental_data = None

    _this_process_label = "XpsWorkChain"

    def update_spectrum_options(self):
        outputs = self._get_child_outputs()
        (
            self.chemical_shifts,
            self.binding_energies,
            self.equivalent_sites_data,
        ) = export_xps_data(outputs)
        options = [key.split("_")[0] for key in self.chemical_shifts.keys()]
        self.spectrum_options = options
        self.spectrum = options[0] if options else None

    def get_data(self):
        if self.spectra_type == "chemical_shift":
            points = self.chemical_shifts
            x_axis_label = "Chemical Shift (eV)"
        else:
            points = self.binding_energies
            x_axis_label = "Binding Energy (eV)"

        self.spectra = xps_spectra_broadening(
            points,
            self.equivalent_sites_data,
            gamma=self.gamma,
            sigma=self.sigma,
            intensity=self.intensity,
        )

        data = [
            {
                "x": d[0],
                "y": d[1],
                "site": site,
            }
            for site, d in self.spectra[self.spectrum].items()
        ]

        fill_type = "tozeroy" if self.fill else None

        return data, x_axis_label, fill_type

    def upload_experimental_data(self, data: dict):
        """Process the uploaded experimental data file."""
        import pandas as pd

        uploaded_file = next(iter(data.values()))
        content = uploaded_file["content"]
        content_str = content.decode("utf-8")

        from io import StringIO

        df = pd.read_csv(StringIO(content_str), header=None)

        self.experimental_data = df
        # Calculate an initial guess for the intensity factor
        total = self.spectra[self._model.spectrum]["total"]
        # Align the max value of the total spectra with the max value of the experimental data
        max_exp = max(self.experimental_data[1])
        max_total = max(total[1])
        self._model.intensity = max_exp / max_total
