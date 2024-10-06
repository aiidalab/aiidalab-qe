import traitlets as tl

from aiida.orm import Group, QueryBuilder, StructureData
from aiidalab_qe.common.panel import SettingsModel

BASE_URL = "https://github.com/superstar54/xps-data/raw/main/pseudo_demo/"


class XpsModel(SettingsModel):
    """Model for the XPS plugin."""

    input_structure = tl.Instance(StructureData, allow_none=True)

    core_hole_treatment = tl.Unicode("xch_smear")
    pseudo_group = tl.Unicode("pseudo_demo_pbe")
    structure_type = tl.Unicode("crystal")
    supercell_min_parameter = tl.Float(8.0)
    calc_binding_energy = tl.Bool(False)
    correction_energies = tl.Dict(
        key_trait=tl.Unicode(),
        value_trait=tl.Dict(),
        default_value={},
    )
    core_levels = tl.Dict(
        key_trait=tl.Unicode(),
        value_trait=tl.Bool(),
        default_value={},
    )

    def update(self):
        if self._pseudo_group_exists():
            self._install_pseudos()

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

        # TODO check logic
        core_level_list = parameters.get("core_level_list", [])
        for orbital in self.core_levels:
            if orbital in core_level_list:
                self.core_levels[orbital] = True  # type: ignore

    def reset(self):
        self.core_hole_treatment = self.traits()["core_hole_treatment"].default_value
        self.pseudo_group = self.traits()["pseudo_group"].default_value
        self.structure_type = self.traits()["structure_type"].default_value
        self.supercell_min_parameter = self.traits()[
            "supercell_min_parameter"
        ].default_value
        self.calc_binding_energy = self.traits()["calc_binding_energy"].default_value

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
