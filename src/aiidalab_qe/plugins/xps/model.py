import traitlets as tl

from aiida.common import NotExistent
from aiida.orm import Group, QueryBuilder, StructureData, load_group
from aiidalab_qe.common.panel import SettingsModel

BASE_URL = "https://github.com/superstar54/xps-data/raw/main/pseudo_demo/"


class XpsModel(SettingsModel):
    dependencies = [
        "input_structure",
    ]

    input_structure = tl.Instance(StructureData, allow_none=True)

    core_hole_treatment = tl.Unicode("xch_smear")
    pseudo_group = tl.Unicode("pseudo_demo_pbe")
    structure_type = tl.Unicode("crystal")
    supercell_min_parameter = tl.Float(8.0)
    calc_binding_energy = tl.Bool(False)
    correction_energies = tl.Dict(
        key_trait=tl.Unicode(),  # <element>_<orbital>
        value_trait=tl.Dict(
            key_trait=tl.Unicode(),
            value_trait=tl.Float(),
        ),
        default_value={},
    )
    core_levels = tl.Dict(
        key_trait=tl.Unicode(),  # core level
        value_trait=tl.Bool(),  # whether the core level is included
        default_value={},
    )

    def update(self, specific=""):
        with self.hold_trait_notifications():
            self._update_correction_energies()
            if not specific or specific == "pseudo_group":
                self._update_pseudos()

    def get_supported_core_levels(self):
        supported_core_levels = {}
        for key in self.correction_energies:
            element = key.split("_")[0]
            if element not in supported_core_levels:
                supported_core_levels[element] = [key]
            else:
                supported_core_levels[element].append(key)
        return supported_core_levels

    def get_model_state(self):
        return {
            # "core_hole_treatment": self.core_hole_treatment,
            "structure_type": self.structure_type,
            "pseudo_group": self.pseudo_group,
            "correction_energies": self.correction_energies,
            "core_level_list": list(self.core_levels.keys()),
        }

    def set_model_state(self, parameters: dict):
        self.pseudo_group = parameters.get(
            "pseudo_group",
            self.traits()["pseudo_group"].default_value,
        )
        self.structure_type = parameters.get(
            "structure_type",
            self.traits()["structure_type"].default_value,
        )

        core_level_list = parameters.get("core_level_list", [])
        for orbital in self.core_levels:
            if orbital in core_level_list:
                self.core_levels[orbital] = True  # type: ignore

    def reset(self):
        with self.hold_trait_notifications():
            for key in [
                "core_hole_treatment",
                "pseudo_group",
                "structure_type",
                "supercell_min_parameter",
                "calc_binding_energy",
            ]:
                setattr(self, key, self.traits()[key].default_value)

    def _update_correction_energies(self):
        try:
            group = load_group(self.pseudo_group)
            self.correction_energies = group.base.extras.get("correction")
        except NotExistent:
            self.correction_energies = {}
            # TODO What if the group does not exist? Should we proceed? Can this happen?

    def _update_pseudos(self):
        if self._pseudo_group_exists():
            self._install_pseudos()

    def _pseudo_group_exists(self, _=None):
        qb = QueryBuilder()
        qb.append(
            Group,
            filters={"label": self.pseudo_group},
        )
        return len(qb.all()) == 0

    def _install_pseudos(self):
        import os
        from pathlib import Path
        from subprocess import run

        url = BASE_URL + self.pseudo_group + ".aiida"

        env = os.environ.copy()
        env["PATH"] = f"{env['PATH']}:{Path.home().joinpath('.local', 'bin')}"

        def run_(*args, **kwargs):
            return run(*args, env=env, capture_output=True, check=True, **kwargs)

        run_(["verdi", "archive", "import", url, "--no-import-group"])
