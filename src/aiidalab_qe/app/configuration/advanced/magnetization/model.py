from copy import deepcopy

import traitlets as tl

from aiidalab_qe.common.mixins import HasInputStructure

from ..subsettings import AdvancedSubModel


class MagnetizationModel(AdvancedSubModel, HasInputStructure):
    dependencies = [
        "input_structure",
        "electronic_type",
        "spin_type",
    ]

    electronic_type = tl.Unicode()
    spin_type = tl.Unicode()
    override = tl.Bool()

    type = tl.Unicode("starting_magnetization")
    total = tl.Float(0.0)
    moments = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
        value_trait=tl.Float(),  # magnetic moment
        default_value={},
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._defaults = {
            "moments": {},
        }

    def update(self, specific=""):  # noqa: ARG002
        if self.spin_type == "none" or self.input_structure is None:
            self._defaults["moments"] = {}
        else:
            self._defaults["moments"] = {
                kind_name: 0.0 for kind_name in self.input_structure.get_kind_names()
            }
        with self.hold_trait_notifications():
            self.moments = self._get_default_moments()

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self.traits()["type"].default_value
            self.total = self.traits()["total"].default_value
            self.moments = self._get_default_moments()

    def _get_default_moments(self):
        return deepcopy(self._defaults["moments"])
