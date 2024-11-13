from aiida import orm
from aiidalab_qe.common.panel import ResultsModel


# TODO if combined, this model should extend `HasModels`, and effectively
# TODO reduce to a container of Bands and PDOS, similar to its results panel
class ElectronicStructureResultsModel(ResultsModel):
    identifier = "electronic_structure"

    def get_pdos_node(self):
        try:
            return self.outputs.pdos
        except AttributeError:
            return None

    def get_bands_node(self):
        # Check for 'bands' or 'bands_projwfc' attributes within 'bands' output
        if self._has_bands:
            return self.outputs.bands.bands
        elif self._has_band_projections:
            return self.outputs.bands.bands_projwfc
        elif self._has_bands_output:
            # If neither 'bands' nor 'bands_projwfc' exist, use 'bands_output' itself
            # This is the case for compatibility with older versions of the plugin
            return self.outputs.bands

    @property
    def include(self):
        # TODO `and` should be `or` if combined
        return all(identifier in self.properties for identifier in ("bands", "pdos"))

    @property
    def has_results(self):
        # TODO first `and` should be `or` if combined
        return self._has_bands and (self._has_dos and self._has_dos_projections)

    @property
    def process_states(self):
        nodes = self._fetch_electronic_structure_process_nodes()
        return [
            node.process_state.value if node and node.process_state else "queued"
            for node in nodes
        ]

    def update_process_state_notification(self):
        processes = ("Bands", "PDOS")
        states = self.process_states
        self.process_state_notification = "\n".join(
            f"""
                <div class="alert alert-{self.CSS_MAP.get(state, "info")}">
                    <b>{process} status:</b> {state.upper()}
                </div>
            """
            for process, state in zip(processes, states)
        )

    @property
    def _has_bands_output(self):
        return hasattr(self.outputs, "bands")

    @property
    def _has_bands(self):
        return self._has_bands_output and hasattr(self.outputs.bands, "bands")

    @property
    def _has_band_projections(self):
        return self._has_bands_output and hasattr(self.outputs.bands, "bands_projwfc")

    @property
    def _has_pdos_output(self):
        return hasattr(self.outputs, "pdos")

    @property
    def _has_dos(self):
        return self._has_pdos_output and hasattr(self.outputs.pdos, "dos")

    @property
    def _has_dos_projections(self):
        return self._has_pdos_output and hasattr(self.outputs.pdos, "projwfc")

    def _fetch_electronic_structure_process_nodes(self):
        return (
            orm.QueryBuilder()
            .append(
                orm.WorkChainNode,
                filters={"uuid": self.process_node.uuid},
                tag="root_process",
            )
            .append(
                orm.WorkChainNode,
                filters={
                    "attributes.process_label": {
                        "in": [
                            "BandsWorkChain",
                            "PdosWorkChain",
                        ]
                    }
                },
                with_incoming="root_process",
            )
            .all(flat=True)
        )
