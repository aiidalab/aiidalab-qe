from copy import deepcopy

import traitlets as tl

from aiida.common.exceptions import NotExistent
from aiida_quantumespresso.workflows.protocols.utils import get_starting_magnetization
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.utils import fetch_pseudo_family_by_label

from ..subsettings import AdvancedCalculationSubSettingsModel


class MagnetizationConfigurationSettingsModel(
    AdvancedCalculationSubSettingsModel,
    HasInputStructure,
):
    identifier = "magnetization"

    dependencies = [
        "input_structure",
        "electronic_type",
        "spin_type",
        "pseudos.family",
    ]

    electronic_type = tl.Unicode()
    spin_type = tl.Unicode()
    family = tl.Unicode()

    type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["Magnetic moments", "starting_magnetization"],
            ["Total magnetization", "tot_magnetization"],
        ],
    )
    type = tl.Unicode("starting_magnetization")
    total = tl.Float(0.0)
    moments = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
        value_trait=tl.Float(),  # magnetic moment
        default_value={},
    )

    def update(self, specific=""):  # noqa: ARG002
        if self.spin_type == "none" or not self.has_structure:
            self._defaults["moments"] = {}
        else:
            try:
                family = fetch_pseudo_family_by_label(self.family)
                initial_guess = get_starting_magnetization(self.input_structure, family)
                self._defaults["moments"] = {
                    kind.name: round(
                        initial_guess[kind.name]
                        * family.get_pseudo(kind.symbol).z_valence,
                        3,
                    )
                    for kind in self.input_structure.kinds
                }
            except NotExistent:
                self._defaults["moments"] = {
                    kind_name: 0.0
                    for kind_name in self.input_structure.get_kind_names()
                }
        with self.hold_trait_notifications():
            self.moments = self._get_default_moments()

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self.traits()["type"].default_value
            self.total = self.traits()["total"].default_value
            self.moments = self._get_default_moments()

    def _get_default_moments(self):
        return deepcopy(self._defaults.get("moments", {}))
