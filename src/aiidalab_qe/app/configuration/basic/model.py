import traitlets as tl

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import PanelModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore

OLD_PROTOCOL_MAP = {
    "moderate": "balanced",
    "precise": "stringent",
}


class BasicConfigurationSettingsModel(
    PanelModel,
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
            ("Balanced", "balanced"),
            ("Stringent", "stringent"),
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

    def get_model_state(self) -> dict:
        return {
            "protocol": self.protocol,
            "spin_type": self.spin_type,
            "electronic_type": self.electronic_type,
        }

    def set_model_state(self, state: dict):
        protocol = state.get("protocol", self.protocol)
        self.protocol = OLD_PROTOCOL_MAP.get(protocol, protocol)
        self.spin_type = state.get("spin_type", self.spin_type)
        self.electronic_type = state.get("electronic_type", self.electronic_type)

    def reset(self):
        with self.hold_trait_notifications():
            self.protocol = self._get_default("protocol")
            self.spin_type = self._get_default("spin_type")
            self.electronic_type = self._get_default("electronic_type")
            self.spin_orbit = self._get_default("spin_orbit")
