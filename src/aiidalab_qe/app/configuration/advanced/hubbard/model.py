from __future__ import annotations

from copy import deepcopy

import numpy as np
import traitlets as tl
from pymatgen.core.periodic_table import Element

from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.common.mixins import HasInputStructure

from ..subsettings import AdvancedCalculationSubSettingsModel


class HubbardConfigurationSettingsModel(
    AdvancedCalculationSubSettingsModel,
    HasInputStructure,
):
    identifier = "hubbard"

    dependencies = [
        "input_structure",
    ]

    is_active = tl.Bool(False)
    has_eigenvalues = tl.Bool(False)
    parameters = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
        value_trait=tl.Float(),  # U value
        default_value={},
    )
    eigenvalue_options = tl.List(
        trait=tl.Unicode(),
        default_value=["-1", "0", "1"],
    )
    eigenvalues = tl.List(
        trait=tl.List(),  # [[[[state, spin, kind name, eigenvalue] # state] # spin] # kind name]
        default_value=[],
    )

    applicable_kind_names = []
    orbital_labels = []

    def update(self, specific=""):  # noqa: ARG002
        if not self.has_structure:
            self.applicable_kind_names = []
            self.orbital_labels = []
            self._defaults |= {
                "parameters": {},
                "eigenvalues": [],
            }
        else:
            self.orbital_labels = self._define_orbital_labels()
            if isinstance(self.input_structure, HubbardStructureData):
                self._defaults["parameters"] = (
                    self.get_parameters_from_hubbard_structure()
                )
                self.is_active = True
            else:
                self._defaults["parameters"] = self._define_default_parameters()
            self.applicable_kind_names = self._define_applicable_kind_names()
            self._defaults["eigenvalues"] = self._define_default_eigenvalues()
        with self.hold_trait_notifications():
            self.parameters = self._get_default_parameters()
            self.eigenvalues = self._get_default_eigenvalues()
            self.needs_eigenvalues_widget = len(self.applicable_kind_names) > 0

    def get_active_eigenvalues(self):
        if not (
            active_eigenvalues := [
                orbital_eigenvalue
                for element_eigenvalues in self.eigenvalues
                for spin_row in element_eigenvalues
                for orbital_eigenvalue in spin_row
                if orbital_eigenvalue[-1] != -1
            ]
        ):
            return []
        eigenvalues_array = np.array(active_eigenvalues, dtype=object)
        new_shape = (int(np.prod(eigenvalues_array.shape[:-1])), 4)
        return [
            tuple(eigenvalue)
            for eigenvalue in eigenvalues_array.reshape(new_shape).tolist()
        ]

    def set_active_eigenvalues(self, eigenvalues: list):
        eigenvalues_array = np.array(eigenvalues, dtype=object)
        num_states = len(set(eigenvalues_array[:, 0]))
        num_spins = len(set(eigenvalues_array[:, 1]))
        num_kinds = len(set(eigenvalues_array[:, 2]))
        new_shape = (num_kinds, num_spins, num_states, 4)
        self.eigenvalues = eigenvalues_array.reshape(new_shape).tolist()
        self.has_eigenvalues = True

    def get_parameters_from_hubbard_structure(self):
        hubbard_parameters = self.input_structure.hubbard.dict()["parameters"]
        sites = self.input_structure.sites
        return {
            f"{sites[hp['atom_index']].kind_name} - {hp['atom_manifold']}": hp["value"]
            for hp in hubbard_parameters
        }

    def reset(self):
        with self.hold_trait_notifications():
            self.is_active = False
            self.has_eigenvalues = False
            self.parameters = self._get_default_parameters()
            self.eigenvalues = self._get_default_eigenvalues()

    def _define_orbital_labels(self):
        hubbard_manifold_list = [
            self._get_manifold(Element(kind.symbol))
            for kind in self.input_structure.kinds
        ]
        return [
            f"{kind_name} - {manifold}"
            for kind_name, manifold in zip(
                self.input_structure.get_kind_names(),
                hubbard_manifold_list,
            )
        ]

    def _define_default_parameters(self):
        return {label: 0.0 for label in self.orbital_labels}

    def _define_applicable_kind_names(self):
        applicable_kind_names = []
        for kind in self.input_structure.kinds:
            element = Element(kind.symbol)
            if (
                element.is_transition_metal
                or element.is_lanthanoid
                or element.is_actinoid
            ):
                num_states = 5 if element.is_transition_metal else 7
                applicable_kind_names.append((kind.name, num_states))
        return applicable_kind_names

    def _define_default_eigenvalues(self):
        return [
            [
                [
                    [state + 1, spin, kind_name, -1]  # default eigenvalue
                    for state in range(num_states)
                ]
                for spin in [1, 2]  # spin up and down
            ]
            for kind_name, num_states in self.applicable_kind_names  # transition metals and lanthanoids
        ]

    def _get_default_parameters(self):
        return deepcopy(self._defaults["parameters"])

    def _get_default_eigenvalues(self):
        return deepcopy(self._defaults["eigenvalues"])

    def _get_manifold(self, element):
        valence = [
            orbital
            for orbital in element.electronic_structure.split(".")
            if "[" not in orbital
        ]
        orbital_shells = [shell[:2] for shell in valence]

        def is_condition_met(shell):
            return condition and condition in shell

        # Conditions for determining the Hubbard manifold
        # to be selected from the electronic structure
        conditions = {
            element.is_transition_metal: "d",
            element.is_lanthanoid or element.is_actinoid: "f",
            element.is_post_transition_metal
            or element.is_metalloid
            or element.is_halogen
            or element.is_chalcogen
            or element.symbol in ["C", "N", "P"]: "p",
            element.is_alkaline or element.is_alkali or element.is_noble_gas: "s",
            element.symbol in ["H", "He"]: "s",
        }

        condition = next(
            (shell for condition, shell in conditions.items() if condition), None
        )

        hubbard_manifold = next(
            (shell for shell in orbital_shells if is_condition_met(shell)), None
        )

        return hubbard_manifold
