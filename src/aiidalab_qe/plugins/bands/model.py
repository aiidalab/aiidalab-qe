import traitlets as tl

from aiidalab_qe.common.panel import PanelModel


class BandsModel(PanelModel):
    """Model for the band structure plugin."""

    kpath_2d = tl.Unicode("hexagonal")

    def get_model_state(self):
        return {
            "kpath_2d": self.kpath_2d,
        }

    def set_model_state(self, parameters: dict):
        self.kpath_2d = parameters.get(
            "kpath_2d",
            self.traits()["kpath_2d"].default_value,
        )

    def reset(self):
        self.kpath_2d = self.traits()["kpath_2d"].default_value
