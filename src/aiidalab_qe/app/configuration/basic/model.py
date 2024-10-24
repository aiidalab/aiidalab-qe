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

    protocol = tl.Unicode(DEFAULT["workchain"]["protocol"])
    spin_type = tl.Unicode(DEFAULT["workchain"]["spin_type"])
    electronic_type = tl.Unicode(DEFAULT["workchain"]["electronic_type"])

    def __init__(self, include=False, *args, **kwargs):
        super().__init__(include, *args, **kwargs)

        self._defaults = {
            "protocol": self.traits()["protocol"].default_value,
            "spin_type": self.traits()["spin_type"].default_value,
            "electronic_type": self.traits()["electronic_type"].default_value,
        }

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
            self.protocol = self._defaults["protocol"]
            self.spin_type = self._defaults["spin_type"]
            self.electronic_type = self._defaults["electronic_type"]
