import traitlets as tl

from aiidalab_qe.common.panel import ConfigurationSettingsModel


class BandsConfigurationSettingsModel(ConfigurationSettingsModel):
    """Model for the band structure plugin."""

    title = "Band structure"
    identifier = "bands"

    projwfc_bands = tl.Bool(False)

    def get_model_state(self):
        return {"projwfc_bands": self.projwfc_bands}

    def set_model_state(self, parameters: dict):
        self.projwfc_bands = parameters.get("projwfc_bands", False)

    def reset(self):
        self.projwfc_bands = False
