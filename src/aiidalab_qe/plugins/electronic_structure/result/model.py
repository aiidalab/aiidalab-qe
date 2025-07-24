from __future__ import annotations

import traitlets as tl

from aiidalab_qe.common.panel import ResultsModel


class ElectronicStructureResultsModel(ResultsModel):
    title = "Electronic Structure"
    identifier = "electronic-structure"

    identifiers = ("bands", "pdos")

    _bands_process_label = "BandsWorkChain"
    _bands_process_uuid = None

    _pdos_process_label = "PdosWorkChain"
    _pdos_process_uuid = None

    _TITLE_MAPPING = {
        "bands": "bands",
        "pdos": "PDOS",
    }

    dos_atoms_group_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("Group by element (atomic species)", "kinds"),
            ("No grouping (each site separately)", "atoms"),
        ],
    )
    dos_atoms_group = tl.Unicode("kinds")
    dos_plot_group_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ("Group all orbitals per atom", "total"),
            ("Group by angular momentum ", "angular_momentum"),
            ("No grouping (each orbital separately)", "orbital"),
        ],
    )
    dos_plot_group = tl.Unicode("angular_momentum")
    selected_atoms = tl.Unicode("")
    project_bands_box = tl.Bool(False)
    proj_bands_width = tl.Float(0.5)

    horizontal_width_percentage = tl.Int(100)
    bands_width_percentage = tl.Int(70)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._completed_processes = set()

    @property
    def include(self):
        return any(identifier in self.properties for identifier in self.identifiers)

    @property
    def has_results(self):
        return self._has_bands or self._has_pdos

    @property
    def needs_property_selector(self):
        return len(self.identifiers) > 1

    def update(self):
        super().update()
        self.identifiers = list(
            filter(
                lambda identifier: identifier in self.properties,
                self.identifiers,
            ),
        )
        parts = [self._TITLE_MAPPING[identifier] for identifier in self.identifiers]
        self.title = f"Electronic {' + '.join(parts)}"

    def update_process_status_notification(self):
        self._status_notifications = []
        for identifier in self.identifiers:
            if identifier in self._completed_processes:
                continue
            status = self._get_child_process_status(identifier)
            self._status_notifications.append(status)
            if "success" in status:
                self._completed_processes.add(identifier)
        self.process_status_notification = "\n".join(self._status_notifications)

    def has_partial_results(self, identifier):
        if identifier == "bands":
            return self._has_bands
        elif identifier == "pdos":
            return self._has_pdos
        return False

    @property
    def _has_bands(self):
        outputs = self._get_child_outputs("bands")
        return any(output in outputs for output in ("bands", "bands_projwfc"))

    @property
    def _has_pdos(self):
        outputs = self._get_child_outputs("pdos")
        return all(output in outputs for output in ("dos", "projwfc"))

    def get_model_state(self):
        return {
            "dos_atoms_group": self.dos_atoms_group,
            "dos_plot_group": self.dos_plot_group,
            "selected_atoms": self.selected_atoms,
            "horizontal_width_percentage": self.horizontal_width_percentage,
            "bands_width_percentage": self.bands_width_percentage,
        }

    def set_model_state(self, parameters):
        self.dos_atoms_group = parameters.get("dos_atoms_group", "kinds")
        self.dos_plot_group = parameters.get("dos_plot_group", "angular_momentum")
        self.selected_atoms = parameters.get("selected_atoms", "")
        self.horizontal_width_percentage = parameters.get(
            "horizontal_width_percentage", 100
        )
        self.bands_width_percentage = parameters.get("bands_width_percentage", 70)
