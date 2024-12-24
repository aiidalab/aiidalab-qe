from __future__ import annotations

import traitlets as tl

from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

from ..subsettings import AdvancedCalculationSubSettingsModel


class SmearingConfigurationSettingsModel(AdvancedCalculationSubSettingsModel):
    identifier = "smearing"

    dependencies = [
        "protocol",
    ]

    protocol = tl.Unicode()

    type_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "cold",
            "gaussian",
            "fermi-dirac",
            "methfessel-paxton",
        ],
    )
    type = tl.Unicode("cold")
    degauss = tl.Float(0.01)

    def update(self, specific=""):  # noqa: ARG002
        parameters = (
            PwBaseWorkChain.get_protocol_inputs(self.protocol)
            .get("pw", {})
            .get("parameters", {})
            .get("SYSTEM", {})
        )
        self._defaults |= {
            "type": parameters["smearing"],
            "degauss": parameters["degauss"],
        }
        with self.hold_trait_notifications():
            self.type = self._defaults["type"]
            self.degauss = self._defaults["degauss"]

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self._get_default("type")
            self.degauss = self._get_default("degauss")

    def _get_default(self, trait):
        return self._defaults.get(trait, self.traits()[trait].default_value)
