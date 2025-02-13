import traitlets as tl

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import ConfigurationSettingsModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class BasicConfigurationSettingsModel(
    ConfigurationSettingsModel,
    HasInputStructure,
):
    title = "Basic settings"
    identifier = "workchain"

    dependencies = [
        "input_structure",
    ]

    protocol_options = tl.List(
        trait=tl.Tuple(tl.Unicode(), tl.Unicode()),
        default_value=[
            ("Fast", "fast"),
            ("Moderate", "moderate"),
            ("Precise", "precise"),
        ],
    )
    protocol = tl.Unicode(DEFAULT["workchain"]["protocol"])
    spin_type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["Off", "none"],
            ["On", "collinear"],
        ],
    )
    spin_type = tl.Unicode(DEFAULT["workchain"]["spin_type"])
    electronic_type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["Metal", "metal"],
            ["Insulator", "insulator"],
        ],
    )
    electronic_type = tl.Unicode(DEFAULT["workchain"]["electronic_type"])
    spin_orbit_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["Off", "wo_soc"],
            ["On", "soc"],
        ],
    )
    spin_orbit = tl.Unicode("wo_soc")

    include = True

    def get_model_state(self):
        return {
            "protocol": self.protocol,
            "spin_type": self.spin_type,
            "electronic_type": self.electronic_type,
        }

    def set_model_state(self, parameters):
        self.protocol = parameters.get("protocol", self.protocol)
        self.spin_type = parameters.get("spin_type", self.spin_type)
        self.electronic_type = parameters.get("electronic_type", self.electronic_type)

    def reset(self):
        with self.hold_trait_notifications():
            self.protocol = self.traits()["protocol"].default_value
            self.spin_type = self.traits()["spin_type"].default_value
            self.electronic_type = self.traits()["electronic_type"].default_value
            self.spin_orbit = self.traits()["spin_orbit"].default_value
