from __future__ import annotations

from copy import deepcopy

import ipywidgets as ipw
import numpy as np
import traitlets as tl
from aiida_pseudo.common.units import U
from pymatgen.core.periodic_table import Element

from aiida import orm
from aiida.common import exceptions
from aiida.plugins import GroupFactory
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.panel import SettingsModel
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION, PseudoFamily

SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class AdvancedModel(SettingsModel):
    dependencies = [
        "input_structure",
        "workchain.protocol",
        "workchain.spin_type",
        "workchain.electronic_type",
    ]

    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    protocol = tl.Unicode()
    spin_type = tl.Unicode()
    electronic_type = tl.Unicode()

    clean_workdir = tl.Bool(False)
    override = tl.Bool(False)
    total_charge = tl.Float(DEFAULT["advanced"]["tot_charge"])
    van_der_waals = tl.Unicode(DEFAULT["advanced"]["vdw_corr"])
    spin_orbit = tl.Unicode("wo_soc")
    forc_conv_thr = tl.Float(0.0)
    forc_conv_thr_step = tl.Float(1e-4)
    etot_conv_thr = tl.Float(0.0)
    etot_conv_thr_step = tl.Float(1e-5)
    scf_conv_thr = tl.Float(0.0)
    scf_conv_thr_step = tl.Float(1e-10)
    electron_maxstep = tl.Int(80)
    kpoints_distance = tl.Float(0.0)
    mesh_grid = tl.Unicode("")

    def __init__(self, include=False, *args, **kwargs):
        super().__init__(include, *args, **kwargs)

        self.dftd3_version = {
            "dft-d3": 3,
            "dft-d3bj": 4,
            "dft-d3m": 5,
            "dft-d3mbj": 6,
        }

        self.smearing = SmearingModel()
        ipw.dlink(
            (self, "protocol"),
            (self.smearing, "protocol"),
        )
        ipw.dlink(
            (self, "override"),
            (self.smearing, "override"),
        )
        self.smearing.observe(
            self._unconfirm,
            tl.All,
        )

        self.magnetization = MagnetizationModel()
        ipw.dlink(
            (self, "input_structure"),
            (self.magnetization, "input_structure"),
        )
        ipw.dlink(
            (self, "electronic_type"),
            (self.magnetization, "electronic_type"),
        )
        ipw.dlink(
            (self, "spin_type"),
            (self.magnetization, "spin_type"),
        )
        ipw.dlink(
            (self, "override"),
            (self.magnetization, "override"),
        )
        self.magnetization.observe(
            self._unconfirm,
            tl.All,
        )

        self.hubbard = HubbardModel()
        ipw.dlink(
            (self, "input_structure"),
            (self.hubbard, "input_structure"),
        )
        ipw.dlink(
            (self, "override"),
            (self.hubbard, "override"),
        )
        self.hubbard.observe(
            self._unconfirm,
            tl.All,
        )

        self.pseudos = PseudosModel()
        ipw.dlink(
            (self, "input_structure"),
            (self.pseudos, "input_structure"),
        )
        ipw.dlink(
            (self, "protocol"),
            (self.pseudos, "protocol"),
        )
        ipw.dlink(
            (self, "spin_orbit"),
            (self.pseudos, "spin_orbit"),
        )
        ipw.dlink(
            (self, "override"),
            (self.pseudos, "override"),
        )
        self.pseudos.observe(
            self._unconfirm,
            tl.All,
        )

        self.observe(
            self._unconfirm,
            tl.All,
        )

        self._defaults = {
            "forc_conv_thr": self.traits()["forc_conv_thr"].default_value,
            "forc_conv_thr_step": self.traits()["forc_conv_thr_step"].default_value,
            "etot_conv_thr": self.traits()["etot_conv_thr"].default_value,
            "etot_conv_thr_step": self.traits()["etot_conv_thr_step"].default_value,
            "scf_conv_thr": self.traits()["scf_conv_thr"].default_value,
            "scf_conv_thr_step": self.traits()["scf_conv_thr_step"].default_value,
            "kpoints_distance": self.traits()["kpoints_distance"].default_value,
        }

    def update(self):
        with self.hold_trait_notifications():
            self._update_defaults()
            self.forc_conv_thr = self._defaults["forc_conv_thr"]
            self.forc_conv_thr_step = self._defaults["forc_conv_thr_step"]
            self.etot_conv_thr = self._defaults["etot_conv_thr"]
            self.etot_conv_thr_step = self._defaults["etot_conv_thr_step"]
            self.scf_conv_thr = self._defaults["scf_conv_thr"]
            self.scf_conv_thr_step = self._defaults["scf_conv_thr_step"]
            self.kpoints_distance = self._defaults["kpoints_distance"]

    def update_kpoints_mesh(self, _=None):
        if self.input_structure is None:
            mesh_grid = ""
        elif self.kpoints_distance > 0:
            mesh = create_kpoints_from_distance.process_class._func(
                self.input_structure,
                orm.Float(self.kpoints_distance),
                orm.Bool(False),
            )
            mesh_grid = f"Mesh {mesh.get_kpoints_mesh()[0]!s}"
        else:
            mesh_grid = "Please select a number higher than 0.0"
        self._defaults["mesh_grid"] = mesh_grid
        self.mesh_grid = mesh_grid

    def get_model_state(self):
        parameters = {
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "tot_charge": self.total_charge,
                    },
                    "CONTROL": {
                        "forc_conv_thr": self.forc_conv_thr,
                        "etot_conv_thr": self.etot_conv_thr,
                    },
                    "ELECTRONS": {
                        "conv_thr": self.scf_conv_thr,
                        "electron_maxstep": self.electron_maxstep,
                    },
                }
            },
            "clean_workdir": self.clean_workdir,
            "pseudo_family": self.pseudos.family,
            "kpoints_distance": self.kpoints_distance,
        }

        if self.hubbard.is_active:
            parameters["hubbard_parameters"] = {"hubbard_u": self.hubbard.parameters}
            if self.hubbard.has_eigenvalues:
                parameters["pw"]["parameters"]["SYSTEM"].update(
                    {"starting_ns_eigenvalue": self.hubbard.get_active_eigenvalues()}
                )

        if self.pseudos.dictionary:
            parameters["pw"]["pseudos"] = self.pseudos.dictionary
            parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = self.pseudos.ecutwfc
            parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"] = self.pseudos.ecutrho

        if self.van_der_waals in ["none", "ts-vdw"]:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = self.van_der_waals
        else:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = "dft-d3"
            parameters["pw"]["parameters"]["SYSTEM"]["dftd3_version"] = (
                self.dftd3_version[self.van_der_waals]
            )

        # there are two choose, use link or parent
        if self.spin_type == "collinear":
            parameters["initial_magnetic_moments"] = self.magnetization.moments
        if self.electronic_type == "metal":
            # smearing type setting
            parameters["pw"]["parameters"]["SYSTEM"]["smearing"] = self.smearing.type
            # smearing degauss setting
            parameters["pw"]["parameters"]["SYSTEM"]["degauss"] = self.smearing.degauss

        # Set tot_magnetization for collinear simulations.
        if self.spin_type == "collinear":
            # Conditions for metallic systems.
            # Select the magnetization type and set the value if override is True
            if self.electronic_type == "metal" and self.override:
                if self.magnetization.type == "tot_magnetization":
                    parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                        self.magnetization.total
                    )
                else:
                    parameters["initial_magnetic_moments"] = self.magnetization.moments
            # Conditions for insulator systems. Default value is 0.0
            elif self.electronic_type == "insulator":
                parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                    self.magnetization.total
                )

        # Spin-Orbit calculation
        if self.spin_orbit == "soc":
            parameters["pw"]["parameters"]["SYSTEM"]["lspinorb"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["noncolin"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["nspin"] = 4

        return parameters

    def set_model_state(self, parameters):
        if "pseudo_family" in parameters:
            pseudo_family = PseudoFamily.from_string(parameters["pseudo_family"])
            library = pseudo_family.library
            accuracy = pseudo_family.accuracy
            self.pseudos.library = f"{library} {accuracy}"
            self.pseudos.functional = pseudo_family.functional

        if "pseudos" in parameters["pw"]:
            self.pseudos.dict = parameters["pw"]["pseudos"]
            self.pseudos.ecutwfc = parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
            self.pseudos.ecutrho = parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"]

        self.kpoints_distance = parameters.get("kpoints_distance", 0.15)

        if (pw_parameters := parameters.get("pw", {}).get("parameters")) is not None:
            self._set_pw_parameters(pw_parameters)

        if magnetic_moments := parameters.get("initial_magnetic_moments"):
            if isinstance(magnetic_moments, (int, float)):
                magnetic_moments = [magnetic_moments]
            if isinstance(magnetic_moments, list):
                magnetic_moments = dict(
                    zip(
                        self.input_structure.get_kind_names(),
                        magnetic_moments,
                    )
                )
            self.magnetization.moments = magnetic_moments

        if parameters.get("hubbard_parameters"):
            self.hubbard.is_active = True
            self.hubbard.parameters = parameters["hubbard_parameters"]["hubbard_u"]
            starting_ns_eigenvalue = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("SYSTEM", {})
                .get("starting_ns_eigenvalue")
            )
            if starting_ns_eigenvalue is not None:
                self.hubbard.has_eigenvalues = True
                self.hubbard.eigenvalues = starting_ns_eigenvalue

    def reset(self):
        with self.hold_trait_notifications():
            self.total_charge = self.traits()["total_charge"].default_value
            self.van_der_waals = self.traits()["van_der_waals"].default_value
            self.forc_conv_thr = self._defaults["forc_conv_thr"]
            self.forc_conv_thr_step = self._defaults["forc_conv_thr_step"]
            self.etot_conv_thr = self._defaults["etot_conv_thr"]
            self.etot_conv_thr_step = self._defaults["etot_conv_thr_step"]
            self.scf_conv_thr = self._defaults["scf_conv_thr"]
            self.scf_conv_thr_step = self._defaults["scf_conv_thr_step"]
            self.electron_maxstep = self.traits()["electron_maxstep"].default_value
            self.spin_orbit = self.traits()["spin_orbit"].default_value
            self.kpoints_distance = self._defaults["kpoints_distance"]
            self.override = self.traits()["override"].default_value

    def _update_defaults(self):
        parameters = PwBaseWorkChain.get_protocol_inputs(self.protocol)

        self._update_kpoints_distance(parameters)
        self.update_kpoints_mesh()

        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        etot_value = num_atoms * parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        self._set_value_and_step("etot_conv_thr", etot_value)

        scf_value = num_atoms * parameters["meta_parameters"]["conv_thr_per_atom"]
        self._set_value_and_step("scf_conv_thr", scf_value)

        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        self._set_value_and_step("forc_conv_thr", forc_value)

    def _update_kpoints_distance(self, parameters):
        kpoints_distance = parameters["kpoints_distance"] if self.has_pbc else 100.0
        self._defaults["kpoints_distance"] = kpoints_distance

    def _set_value_and_step(self, attribute, value):
        self._defaults[attribute] = value
        if value != 0:
            order_of_magnitude = np.floor(np.log10(abs(value)))
            step = 10 ** (order_of_magnitude - 1)
        else:
            step = 0.1
        self._defaults[f"{attribute}_step"] = step

    def _set_pw_parameters(self, pw_parameters):
        system_params = pw_parameters.get("SYSTEM", {})
        control_params = pw_parameters.get("CONTROL", {})
        electron_params = pw_parameters.get("ELECTRONS", {})

        self.forc_conv_thr = control_params.get("forc_conv_thr", 0.0)
        self.etot_conv_thr = control_params.get("etot_conv_thr", 0.0)
        self.scf_conv_thr = electron_params.get("conv_thr", 0.0)
        self.electron_maxstep = electron_params.get("electron_maxstep", 80)

        self.total_charge = system_params.get("tot_charge", 0)
        self.spin_orbit = "soc" if "lspinorb" in system_params else "wo_soc"

        self.van_der_waals = self.dftd3_version.get(
            system_params.get("dftd3_version"),
            system_params.get("vdw_corr", "none"),
        )

        if "degauss" in system_params:
            self.smearing.degauss = system_params["degauss"]

        if "smearing" in system_params:
            self.smearing.type = system_params["smearing"]

        if "tot_magnetization" in system_params:
            self.magnetization.type = "tot_magnetization"


class AdvancedSubModel(tl.HasTraits):
    _defaults = {}

    def update(self):
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

    def _update_defaults(self):
        raise NotImplementedError


class SmearingModel(AdvancedSubModel):
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

    def update(self):
        with self.hold_trait_notifications():
            self._update_defaults()
            self.type = self._defaults["type"]
            self.degauss = self._defaults["degauss"]

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self._defaults["type"]
            self.degauss = self._defaults["degauss"]

    def _update_defaults(self):
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


class MagnetizationModel(AdvancedSubModel):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    electronic_type = tl.Unicode()
    spin_type = tl.Unicode()
    override = tl.Bool()

    type = tl.Unicode("starting_magnetization")
    total = tl.Float(0.0)
    moments = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # magnetic moment
        default_value={},
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._defaults = {
            "moments": {},
        }

    def update(self):
        with self.hold_trait_notifications():
            self._update_defaults()
            self.moments = self._get_default_moments()

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self.traits()["type"].default_value
            self.total = self.traits()["total"].default_value
            self.moments = self._get_default_moments()

    def _update_defaults(self):
        if self.spin_type == "none" or self.input_structure is None:
            self._defaults["moments"] = {}
        else:
            self._defaults["moments"] = {
                symbol: 0.0 for symbol in self.input_structure.get_kind_names()
            }

    def _get_default_moments(self):
        return deepcopy(self._defaults["moments"])


class HubbardModel(AdvancedSubModel):
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

    applicable_elements = []
    orbital_labels = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._defaults = {
            "parameters": {},
            "eigenvalues": [],
        }

    def update(self):
        with self.hold_trait_notifications():
            self._update_defaults()
            self.parameters = self._get_default_parameters()
            self.eigenvalues = self._get_default_eigenvalues()
            self.needs_eigenvalues_widget = len(self.applicable_elements) > 0

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

    def _update_defaults(self):
        if self.input_structure is None:
            self.applicable_elements = []
            self.orbital_labels = []
            self._defaults.update(
                {
                    "parameters": {},
                    "eigenvalues": [],
                }
            )
        else:
            self.orbital_labels = self._get_labels()
            self._defaults["parameters"] = {label: 0.0 for label in self.orbital_labels}
            self.applicable_elements = [
                *filter(
                    lambda element: (
                        element.is_transition_metal
                        or element.is_lanthanoid
                        or element.is_actinoid
                    ),
                    [
                        Element(symbol)
                        for symbol in self.input_structure.get_kind_names()
                    ],
                )
            ]
            self._defaults["eigenvalues"] = [
                [
                    [
                        [state + 1, spin, element.symbol, -1]  # default eigenvalue
                        for state in range(5 if element.is_transition_metal else 7)
                    ]
                    for spin in range(2)  # spin up and down
                ]
                for element in self.applicable_elements  # transition metals and lanthanoids
            ]

    def _get_default_parameters(self):
        return deepcopy(self._defaults["parameters"])

    def _get_default_eigenvalues(self):
        return deepcopy(self._defaults["eigenvalues"])

    def _get_labels(self):
        symbols = self.input_structure.get_kind_names()
        hubbard_manifold_list = [
            self._get_manifold(Element(symbol))
            for symbol in self.input_structure.get_kind_names()
        ]
        return [
            f"{symbol} - {manifold}"
            for symbol, manifold in zip(symbols, hubbard_manifold_list)
        ]

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


class PseudosModel(AdvancedSubModel):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    protocol = tl.Unicode()
    spin_orbit = tl.Unicode()
    override = tl.Bool()

    dictionary = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Unicode(),  # pseudopotential node uuid
        default_value={},
    )
    family = tl.Unicode(
        "/".join(
            [
                DEFAULT["advanced"]["pseudo_family"]["library"],
                str(DEFAULT["advanced"]["pseudo_family"]["version"]),
                DEFAULT["advanced"]["pseudo_family"]["functional"],
                DEFAULT["advanced"]["pseudo_family"]["accuracy"],
            ]
        )
    )
    functional = tl.Unicode(DEFAULT["advanced"]["pseudo_family"]["functional"])
    functional_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "PBE",
            "PBEsol",
        ],
    )
    library = tl.Unicode(
        " ".join(
            [
                DEFAULT["advanced"]["pseudo_family"]["library"],
                DEFAULT["advanced"]["pseudo_family"]["accuracy"],
            ]
        )
    )
    library_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "SSSP efficiency",
            "SSSP precision",
            "PseudoDojo standard",
            "PseudoDojo stringent",
        ],
    )
    cutoffs = tl.List(
        trait=tl.List(tl.Float()),  # [[ecutwfc values], [ecutrho values]]
        default_value=[[0.0], [0.0]],
    )
    ecutwfc = tl.Float()
    ecutrho = tl.Float()
    status_message = tl.Unicode()
    family_help_message = tl.Unicode()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        ipw.dlink(
            (self, "cutoffs"),
            (self, "ecutwfc"),
            lambda cutoffs: max(cutoffs[0]),
        )
        ipw.dlink(
            (self, "cutoffs"),
            (self, "ecutrho"),
            lambda cutoffs: max(cutoffs[1]),
        )

        self.PSEUDO_HELP_SOC = """
            <div class="pseudo-text">
                Spin-orbit coupling (SOC) calculations are supported exclusively with
                PseudoDojo pseudopotentials. PseudoDojo offers these pseudopotentials
                in two versions: standard and stringent. Here, we utilize the FR
                (fully relativistic) type from PseudoDojo. Please ensure you choose
                appropriate cutoff values for your calculations.
            </div>
        """

        self.PSEUDO_HELP_WO_SOC = """
            <div class="pseudo-text">
                If you are unsure, select 'SSSP efficiency', which for most
                calculations will produce sufficiently accurate results at
                comparatively small computational costs. If your calculations require a
                higher accuracy, select 'SSSP accuracy' or 'PseudoDojo stringent',
                which will be computationally more expensive. SSSP is the standard
                solid-state pseudopotentials. The PseudoDojo used here has the SR
                relativistic type.
            </div>
        """

        self.family_help_message = self.PSEUDO_HELP_WO_SOC

        self._defaults = {
            "family": self.traits()["family"].default_value,
            "functional": self.traits()["functional"].default_value,
            "functional_options": [
                "PBE",
                "PBEsol",
            ],
            "library": self.traits()["library"].default_value,
            "library_options": [
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ],
            "dictionary": {},
            "cutoffs": [[0.0], [0.0]],
        }

    def update(self):
        with self.hold_trait_notifications():
            self._update_defaults()

    def update_default_pseudos(self):
        try:
            pseudo_family = self._get_pseudo_family_from_database()
            pseudos = pseudo_family.get_pseudos(structure=self.input_structure)
        except ValueError as exception:
            self._status_message = f"""
                <div class='alert alert-danger'>
                    ERROR: {exception!s}
                </div>
            """
            return

        self._defaults["dictionary"] = {
            kind: pseudo.uuid for kind, pseudo in pseudos.items()
        }
        self.dictionary = self._get_default_dictionary()

    def update_default_cutoffs(self):
        """Update wavefunction and density cutoffs from pseudo family."""
        try:
            pseudo_family = self._get_pseudo_family_from_database()
            current_unit = pseudo_family.get_cutoffs_unit()
            cutoff_dict = pseudo_family.get_cutoffs()
        except exceptions.NotExistent:
            self._status_message = f"""
                <div class='alert alert-danger'>
                    ERROR: required pseudo family `{self.family}` is
                    not installed. Please use `aiida-pseudo install` to install
                    it."
                </div>
            """
        except ValueError as exception:
            self._status_message = f"""
                <div class='alert alert-danger'>
                    ERROR: failed to obtain recommended cutoffs for pseudos
                    `{pseudo_family}`: {exception}
                </div>
            """
        else:
            symbols = (
                self.input_structure.get_kind_names() if self.input_structure else []
            )

        ecutwfc_list = []
        ecutrho_list = []
        for symbol in symbols:
            cutoff = cutoff_dict.get(symbol, {})
            ecutrho, ecutwfc = (
                U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                for v in cutoff.values()
            )
            ecutwfc_list.append(ecutwfc)
            ecutrho_list.append(ecutrho)

        self._defaults["cutoffs"] = [ecutwfc_list or [0.0], ecutrho_list or [0.0]]
        self.cutoffs = self._get_default_cutoffs()

    def update_library_options(self):
        if self.spin_orbit == "soc":
            library_options = [
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.family_help_message = self.PSEUDO_HELP_SOC
        else:
            library_options = [
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.family_help_message = self.PSEUDO_HELP_WO_SOC
        self._defaults["library_options"] = library_options
        self.library_options = self._defaults["library_options"]

        self.update_family_parameters()

    def update_family_parameters(self):
        if self.spin_orbit == "soc":
            if self.protocol in ["fast", "moderate"]:
                pseudo_family_string = "PseudoDojo/0.4/PBE/FR/standard/upf"
            else:
                pseudo_family_string = "PseudoDojo/0.4/PBE/FR/stringent/upf"
        else:
            pseudo_family_string = PwBaseWorkChain.get_protocol_inputs(self.protocol)[
                "pseudo_family"
            ]

        pseudo_family = PseudoFamily.from_string(pseudo_family_string)

        self._defaults["library"] = f"{pseudo_family.library} {pseudo_family.accuracy}"
        self._defaults["functional"] = pseudo_family.functional

        with self.hold_trait_notifications():
            self.library = self._defaults["library"]
            self.functional = self._defaults["functional"]

    def update_family(self):
        library, accuracy = self.library.split()
        functional = self.functional
        # XXX (jusong.yu): a validator is needed to check the family string is
        # consistent with the list of pseudo families defined in the setup_pseudos.py
        if library == "PseudoDojo":
            if self.spin_orbit == "soc":
                pseudo_family_string = (
                    f"PseudoDojo/{PSEUDODOJO_VERSION}/{functional}/FR/{accuracy}/upf"
                )
            else:
                pseudo_family_string = (
                    f"PseudoDojo/{PSEUDODOJO_VERSION}/{functional}/SR/{accuracy}/upf"
                )
        elif library == "SSSP":
            pseudo_family_string = f"SSSP/{SSSP_VERSION}/{functional}/{accuracy}"
        else:
            raise ValueError(
                f"Unknown pseudo family parameters: {library} | {accuracy}"
            )

        self._defaults["family"] = pseudo_family_string
        self.family = self._defaults["family"]

    def reset(self):
        with self.hold_trait_notifications():
            self.dictionary = self._get_default_dictionary()
            self.cutoffs = self._get_default_cutoffs()
            self.library_options = self._defaults["library_options"]
            self.library = self._defaults["library"]
            self.functional = self._defaults["functional"]
            self.functional_options = self._defaults["functional_options"]
            self.family = self._defaults["family"]
            self.family_help_message = self.PSEUDO_HELP_WO_SOC
            self.status_message = ""

    def _update_defaults(self):
        if self.input_structure is None:
            self._defaults.update(
                {
                    "dictionary": {},
                    "cutoffs": [[0.0], [0.0]],
                }
            )
        else:
            self.update_default_pseudos()
            self.update_default_cutoffs()
        self.update_family_parameters()
        self.update_family()

    def _get_pseudo_family_from_database(self):
        """Get the pseudo family from the database."""
        return (
            orm.QueryBuilder()
            .append(
                (
                    PseudoDojoFamily,
                    SsspFamily,
                    CutoffsPseudoPotentialFamily,
                ),
                filters={"label": self.family},
            )
            .one()[0]
        )

    def _get_default_dictionary(self):
        return deepcopy(self._defaults["dictionary"])

    def _get_default_cutoffs(self):
        return deepcopy(self._defaults["cutoffs"])
