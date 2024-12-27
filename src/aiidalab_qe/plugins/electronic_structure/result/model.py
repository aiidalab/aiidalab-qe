from __future__ import annotations

from aiidalab_qe.common.panel import ResultsModel


# TODO if combined, this model should extend `HasModels`, and effectively
# TODO reduce to a container of Bands and PDOS, similar to its results panel
class ElectronicStructureResultsModel(ResultsModel):
    title = "Electronic Structure"
    identifier = "electronic_structure"

    identifiers = ("bands", "pdos")

    _bands_process_label = "BandsWorkChain"
    _bands_process_uuid = None

    _pdos_process_label = "PdosWorkChain"
    _pdos_process_uuid = None

    _TITLE_MAPPING = {
        "bands": "bands",
        "pdos": "PDOS",
    }

    @property
    def include(self):
        return all(identifier in self.properties for identifier in self.identifiers)

    @property
    def has_results(self):
        return self._has_bands and self._has_pdos

    def update(self):
        super().update()
        electronic_properties = [
            self._TITLE_MAPPING[identifier]
            for identifier in self.identifiers
            if identifier in self.properties
        ]
        self.title = f"Electronic {' + '.join(electronic_properties)}"

    def update_process_status_notification(self):
        statuses = [self._get_child_process_status(child) for child in self.identifiers]
        self.process_status_notification = "\n".join(statuses)

    @property
    def _has_bands(self):
        outputs = self._get_child_outputs("bands")
        return any(output in outputs for output in ("bands", "bands_projwfc"))

    @property
    def _has_pdos(self):
        outputs = self._get_child_outputs("pdos")
        return all(output in outputs for output in ("dos", "projwfc"))
