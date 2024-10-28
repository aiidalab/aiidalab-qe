import traitlets as tl

from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

from ..subsettings import AdvancedSubModel


class SmearingModel(AdvancedSubModel):
    dependencies = [
        "protocol",
    ]

    protocol = tl.Unicode()
    override = tl.Bool()

    type = tl.Unicode("cold")
    degauss = tl.Float(0.0)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._defaults = {
            "type": self.traits()["type"].default_value,
            "degauss": self.traits()["degauss"].default_value,
        }

    def update(self, which):
        with self.hold_trait_notifications():
            self._update_defaults(which)
            self.type = self._defaults["type"]
            self.degauss = self._defaults["degauss"]

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self._defaults["type"]
            self.degauss = self._defaults["degauss"]

    def _update_defaults(self, which):
        parameters = (
            PwBaseWorkChain.get_protocol_inputs(self.protocol)
            .get("pw", {})
            .get("parameters", {})
            .get("SYSTEM", {})
        )
        self._defaults.update(
            {
                "type": parameters["smearing"],
                "degauss": parameters["degauss"],
            }
        )
