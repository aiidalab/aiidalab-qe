from __future__ import annotations

from copy import deepcopy

import traitlets as tl
from pymatgen.core.periodic_table import Element

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

from ..subsettings import AdvancedSubModel


class HubbardModel(AdvancedSubModel):
    dependencies = [
        "input_structure",
    ]

    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    override = tl.Bool()

    is_active = tl.Bool(False)
    has_eigenvalues = tl.Bool(False)
    parameters = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # U value
        default_value={},
    )
    eigenvalues = tl.List(
        trait=tl.List(),  # [[[[state, spin, symbol, eigenvalue] # state] # spin] # symbol]
        default_value=[],
    )

    applicable_kinds = []
    orbital_labels = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._defaults = {
            "parameters": {},
            "eigenvalues": [],
        }

    def update(self, which):
        with self.hold_trait_notifications():
            self._update_defaults(which)
            self.parameters = self._get_default_parameters()
            self.eigenvalues = self._get_default_eigenvalues()
            self.needs_eigenvalues_widget = len(self.applicable_kinds) > 0

    def get_active_eigenvalues(self):
        return [
            orbital_eigenvalue
            for element_eigenvalues in self.eigenvalues
            for spin_row in element_eigenvalues
            for orbital_eigenvalue in spin_row
            if orbital_eigenvalue[-1] != -1
        ]

    def set_parameters_from_hubbard_structure(self):
        hubbard_parameters = self.input_structure.hubbard.dict()["parameters"]
        sites = self.input_structure.sites
        parameters = {
            f"{sites[hp['atom_index']].kind_name} - {hp['atom_manifold']}": hp["value"]
            for hp in hubbard_parameters
        }
        with self.hold_trait_notifications():
            self.parameters = parameters
            self.is_active = True

    def reset(self):
        with self.hold_trait_notifications():
            self.is_active = False
            self.has_eigenvalues = False
            self.parameters = self._get_default_parameters()
            self.eigenvalues = self._get_default_eigenvalues()

    def _update_defaults(self, which):
        if self.input_structure is None:
            self.applicable_kinds = []
            self.orbital_labels = []
            self._defaults |= {
                "parameters": {},
                "eigenvalues": [],
            }
        else:
            self.orbital_labels = self._define_orbital_labels()
            self._defaults["parameters"] = self._define_default_parameters()
            self.applicable_kinds = self._define_applicable_kinds()
            self._defaults["eigenvalues"] = self._define_default_eigenvalues()

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

    def _define_applicable_kinds(self):
        applicable_kinds = []
        for kind in self.input_structure.kinds:
            element = Element(kind.symbol)
            if (
                element.is_transition_metal
                or element.is_lanthanoid
                or element.is_actinoid
            ):
                num_states = 5 if element.is_transition_metal else 7
                applicable_kinds.append((kind, num_states))
        return applicable_kinds

    def _define_default_eigenvalues(self):
        return [
            [
                [
                    [state + 1, spin, kind.symbol, -1]  # default eigenvalue
                    for state in range(num_states)
                ]
                for spin in range(2)  # spin up and down
            ]
            for kind, num_states in self.applicable_kinds  # transition metals and lanthanoids
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
        }

        condition = next(
            (shell for condition, shell in conditions.items() if condition), None
        )

        hubbard_manifold = next(
            (shell for shell in orbital_shells if is_condition_met(shell)), None
        )

        return hubbard_manifold
