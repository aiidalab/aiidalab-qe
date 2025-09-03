import traitlets as tl

from aiidalab_qe.common.panel import PanelModel


class BandsConfigurationSettingsModel(PanelModel):
    """Model for the band structure plugin."""

    title = "Band structure"
    identifier = "bands"

    projwfc_bands = tl.Bool(False)

    def get_model_state(self) -> dict:
        return {"projwfc_bands": self.projwfc_bands}

    def set_model_state(self, state: dict):
        self.projwfc_bands = state.get("projwfc_bands", False)

    def reset(self):
        self.projwfc_bands = False
