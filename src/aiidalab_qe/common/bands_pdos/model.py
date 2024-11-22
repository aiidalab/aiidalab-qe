import base64
import json

import ipywidgets as ipw
import numpy as np
import traitlets as tl
from IPython.display import display

from aiida.common.extendeddicts import AttributeDict
from aiidalab_qe.common.bands_pdos.utils import (
    get_bands_data,
    get_bands_projections_data,
    get_pdos_data,
)
from aiidalab_qe.common.mvc import Model
from aiidalab_widgets_base.utils import string_range_to_list

from .bandpdosplotly import BandsPdosPlotly


class BandsPdosModel(Model):
    bands = tl.Instance(AttributeDict, allow_none=True)
    pdos = tl.Instance(AttributeDict, allow_none=True)

    dos_atoms_group_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("Kinds", "kinds"),
            ("Atomic position", "atoms"),
        ],
    )
    dos_atoms_group = tl.Unicode("kinds")
    dos_plot_group_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("Total", "total"),
            ("Orbital", "orbital"),
            ("Angular momentum", "angular_momentum"),
        ],
    )
    dos_plot_group = tl.Unicode("total")
    selected_atoms = tl.Unicode("")
    project_bands_box = tl.Bool(False)
    proj_bands_width = tl.Float(0.5)

    needs_pdos_options = tl.Bool(False)
    needs_projections_controls = tl.Bool(False)

    pdos_data = {}
    bands_data = {}
    bands_projections_data = {}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        ipw.dlink(
            (self, "bands"),
            (self, "needs_projections_controls"),
            lambda _: self._has_bands_projections,
        )
        ipw.dlink(
            (self, "pdos"),
            (self, "needs_pdos_options"),
            lambda _: self._has_pdos or self.needs_projections_controls,
        )
        ipw.dlink(
            (self, "needs_projections_controls"),
            (self, "needs_pdos_options"),
            lambda _: self._has_pdos or self.needs_projections_controls,
        )

    def fetch_data(self):
        """Fetch the data from the nodes."""
        if self.bands:
            if not self.bands_data:
                self.bands_data = self._get_bands_data()

        if self.pdos:
            self.pdos_data = self._get_pdos_data()

    def create_plot(self):
        """Create the plot."""
        self.helper = BandsPdosPlotly(
            bands_data=self.bands_data,
            bands_projections_data=None,
            pdos_data=self.pdos_data,
        )
        self.plot = self.helper.bandspdosfigure

    def update_bands_projections(self):
        """Update the bands projections."""
        if self.project_bands_box:
            if self.bands_projections_data:
                self._remove_bands_traces()
            self.bands_projections_data = self._get_bands_projections_data()
            self.helper.project_bands = self.bands_projections_data
            self.helper.adding_projected_bands(self.plot)
        else:
            self._remove_bands_traces()

    def _remove_bands_traces(self):
        """Remove the bands traces."""
        self.plot.data = tuple(
            trace
            for trace in self.plot.data
            if trace.xaxis == "x2" or trace.legendgroup == ""
        )
        self.bands_projections_data = {}
        self.helper.project_bands = {}

    def _remove_pdos_traces(self):
        """Remove the PDOS traces."""
        if self.helper.plot_type == "pdos":
            self.plot.data = ()
        elif self.helper.plot_type == "combined":
            self.plot.data = tuple(
                trace for trace in self.plot.data if trace.xaxis != "x2"
            )

    def update_pdos_plot(self):
        """Update the PDOS plot."""
        self.fetch_data()
        if self.project_bands_box:
            self.update_bands_projections()
        self.helper.pdos_data = self.pdos_data
        self._remove_pdos_traces()
        self.helper.adding_pdos_traces(self.plot)

    @property
    def _has_bands(self):
        return bool(self.bands)

    @property
    def _has_pdos(self):
        return bool(self.pdos)

    @property
    def _has_bands_projections(self):
        return self._has_bands and "projwfc" in self.bands

    def _get_pdos_data(self):
        if not self.pdos:
            return None
        expanded_selection, syntax_ok = string_range_to_list(
            self.selected_atoms, shift=-1
        )
        if syntax_ok:
            return get_pdos_data(
                self.pdos,
                group_tag=self.dos_atoms_group,
                plot_tag=self.dos_plot_group,
                selected_atoms=expanded_selection,
            )
        return None

    def _get_bands_data(self):
        if not self.bands:
            return None

        bands_data = get_bands_data(self.bands)
        return bands_data

    def _get_bands_projections_data(self):
        if not self.bands:
            return None

        expanded_selection, syntax_ok = string_range_to_list(
            self.selected_atoms, shift=-1
        )
        if syntax_ok:
            bands_projections = get_bands_projections_data(
                self.bands,
                bands_data=self.bands_data,
                group_tag=self.dos_atoms_group,
                plot_tag=self.dos_plot_group,
                selected_atoms=expanded_selection,
                bands_width=self.proj_bands_width,
            )
            return bands_projections
        return None

    def download_data(self, _=None):
        """Function to download the data."""
        if self.bands_data:
            bands_data_export = {
                key: value.tolist() if isinstance(value, np.ndarray) else value
                for key, value in self.bands_data.items()
            }
            json_str = json.dumps(bands_data_export)
            b64_str = base64.b64encode(json_str.encode()).decode()
            file_name_bands = "bands_data.json"
            self._download(payload=b64_str, filename=file_name_bands)
        if self.pdos_data:
            json_str = json.dumps(self.pdos_data)
            b64_str = base64.b64encode(json_str.encode()).decode()
            file_name_pdos = "dos_data.json"
            self._download(payload=b64_str, filename=file_name_pdos)

    @staticmethod
    def _download(payload, filename):
        """Download payload as a file named as filename."""
        from IPython.display import Javascript

        javas = Javascript(
            f"""
            var link = document.createElement('a');
            link.href = 'data:text/json;charset=utf-8;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """
        )
        display(javas)
