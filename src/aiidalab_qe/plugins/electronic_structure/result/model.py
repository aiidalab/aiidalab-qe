from __future__ import annotations

from aiida.common.extendeddicts import AttributeDict
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
        if not (node := self._fetch_child_process_node("pdos")):
            return
        outputs = {key: getattr(node.outputs, key) for key in node.outputs}
        return AttributeDict(outputs)

    def get_bands_node(self):
        if not (node := self._fetch_child_process_node("bands")):
            return
        if self._has_bands:
            return node.outputs.bands
        elif self._has_band_projections:
            return node.outputs.bands_projwfc
        else:
            # If neither 'bands' nor 'bands_projwfc' exist, use the 'outputs' node itself
            # This is the case for compatibility with older versions of the plugin
            return node.outputs

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
        bands_node = self._fetch_child_process_node("bands")
        pdos_node = self._fetch_child_process_node("pdos")
        return [
            node.process_state.value if node and node.process_state else "queued"
            for node in (bands_node, pdos_node)
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
        node = self._fetch_child_process_node("bands")
        return node and "bands" in node.outputs

    @property
    def _has_band_projections(self):
        node = self._fetch_child_process_node("bands")
        return node and "bands_projwfc" in node.outputs

    @property
    def _has_dos(self):
        node = self._fetch_child_process_node("pdos")
        return node and "dos" in node.outputs

    @property
    def _has_dos_projections(self):
        node = self._fetch_child_process_node("pdos")
        return node and "projwfc" in node.outputs
