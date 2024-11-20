import traitlets as tl

from aiida import orm
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.panel import SettingsModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class WorkChainModel(SettingsModel):
    dependencies = [
        "input_structure",
    ]

    input_structure = tl.Union([tl.Instance(orm.StructureData)], allow_none=True)

    protocol_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "fast",
            "moderate",
            "precise",
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
