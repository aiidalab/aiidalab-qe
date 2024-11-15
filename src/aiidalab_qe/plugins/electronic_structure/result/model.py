from __future__ import annotations

from aiidalab_qe.common.panel import ResultsModel


# TODO if combined, this model should extend `HasModels`, and effectively
# TODO reduce to a container of Bands and PDOS, similar to its results panel
class ElectronicStructureResultsModel(ResultsModel):
    identifier = "electronic_structure"

    _bands_process_label = "BandsWorkChain"
    _bands_process_uuid = None

    _pdos_process_label = "PdosWorkChain"
    _pdos_process_uuid = None

    def get_pdos_node(self):
        return self._get_child_outputs("pdos")

    def get_bands_node(self):
        outputs = self._get_child_outputs("bands")
        if "bands" in outputs:
            return outputs.bands
        elif "bands_projwfc" in outputs:
            return outputs.bands_projwfc
        else:
            # If neither 'bands' nor 'bands_projwfc' exist, use 'bands_output' itself
            # This is the case for compatibility with older versions of the plugin
            return outputs

    @property
    def include(self):
        return all(identifier in self.properties for identifier in ("bands", "pdos"))

    @property
    def has_results(self):
        # TODO should be `or` if merged
        return self._has_bands and self._has_pdos

    @property
    def process_states(self):
        return [
            node.process_state.value if node and node.process_state else "queued"
            for node in (
                self._fetch_child_process_node("bands"),
                self._fetch_child_process_node("pdos"),
            )
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
    def _has_bands(self):
        outputs = self._get_child_outputs("bands")
        return any(output in outputs for output in ("bands", "bands_projwfc"))

    @property
    def _has_pdos(self):
        outputs = self._get_child_outputs("pdos")
        return all(output in outputs for output in ("dos", "projwfc"))
