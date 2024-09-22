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
from aiida_quantumespresso.common.types import RelaxType
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.panel import PanelModel
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION, PseudoFamily

SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class WorkChainModel(PanelModel):
    protocol = tl.Unicode(DEFAULT["workchain"]["protocol"])
    relax_type = tl.Unicode("positions_cell")
    spin_type = tl.Unicode(DEFAULT["workchain"]["spin_type"])
    electronic_type = tl.Unicode(DEFAULT["workchain"]["electronic_type"])

    def get_model_state(self):
        return {
            "protocol": self.protocol,
            "relax_type": self.relax_type,
            "spin_type": self.spin_type,
            "electronic_type": self.electronic_type,
        }

    def set_model_state(self, parameters):
        """Update the settings based on the given dict."""
        for key in [
            "relax_type",
            "spin_type",
            "electronic_type",
        ]:
            if key in parameters:
                setattr(self, key, parameters[key])
        if "protocol" in parameters:
            self.protocol = parameters["protocol"]

    def reset(self):
        with self.hold_trait_notifications():
            self.protocol = self.traits()["protocol"].default_value
            self.relax_type = self.traits()["relax_type"].default_value
            self.spin_type = self.traits()["spin_type"].default_value
            self.electronic_type = self.traits()["electronic_type"].default_value


class SmearingModel(tl.HasTraits):
    protocol = tl.Unicode()
    override = tl.Bool()

    type = tl.Unicode()
    degauss = tl.Float()

    _defaults = {}

    def set_defaults_from_protocol(self):
        parameters = (
            PwBaseWorkChain.get_protocol_inputs(self.protocol)
            .get("pw", {})
            .get("parameters", {})
            .get("SYSTEM", {})
        )
        self._defaults = {
            "type": parameters["smearing"],
            "degauss": parameters["degauss"],
        }

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self._get_default_type()
            self.degauss = self._get_default_degauss()

    @tl.default("type")
    def _get_default_type(self):
        return self._defaults["type"]

    @tl.default("degauss")
    def _get_default_degauss(self):
        return self._defaults["degauss"]


class MagnetizationModel(tl.HasTraits):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    electronic_type = tl.Unicode()
    override = tl.Bool()

    type = tl.Unicode("starting_magnetization")
    total = tl.Float(0.0)
    moments = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # magnetic moment
        default_value={},
    )

    _default_moments = {}

    def set_defaults_from_structure(self):
        if self.input_structure is None:
            self._default_moments = {}
        else:
            self._default_moments = {
                kind.symbol: 0.0 for kind in self.input_structure.kinds
            }
        self.moments = self._get_moments()

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self.traits()["type"].default_value
            self.total = self.traits()["total"].default_value
            self.moments = self._get_moments()

    def _get_moments(self):
        return deepcopy(self._default_moments)


class HubbardModel(tl.HasTraits):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    override = tl.Bool()

    activate = tl.Bool(False)
    eigenvalues_label = tl.Bool(False)  # TODO: rename (widget also)
    parameters = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # U value
        default_value={},
    )
    eigenvalues = tl.List(
        trait=tl.List(),  # [[[[state, spin, kind, eigenvalue] # state] # spin] # kind]
        default_value=[],
    )

    _default_parameters = {}
    _default_eigenvalues = []

    def set_defaults_from_structure(self):
        if self.input_structure is None:
            self.elements = []  # TODO: rename for clarity
            self.input_labels = []  # TODO: rename for clarity
            self._default_parameters = {}
            self._default_eigenvalues = []
        else:
            self.input_labels = self._get_labels()
            self._default_parameters = {label: 0.0 for label in self.input_labels}
            self.elements = [
                *filter(
                    lambda element: (
                        element.is_transition_metal
                        or element.is_lanthanoid
                        or element.is_actinoid
                    ),
                    [Element(kind.symbol) for kind in self.input_structure.kinds],
                )
            ]
            self._default_eigenvalues = [
                [
                    [
                        [state + 1, spin, element.symbol, "-1"]  # default eigenvalue
                        for state in range(5 if element.is_transition_metal else 7)
                    ]
                    for spin in range(2)  # spin up and down
                ]
                for element in self.elements  # transition metals and lanthanoids
            ]
        self.parameters = self._get_default_parameters()
        self.eigenvalues = self._get_default_eigenvalues()
        self.needs_eigenvalues_widget = len(self.elements) > 0

    def set_parameters_from_hubbard_structure(self):
        hubbard_parameters = self.input_structure.hubbard.dict()["parameters"]
        sites = self.input_structure.sites
        parameters = {
            f"{sites[hp['atom_index']].kind_name} - {hp['atom_manifold']}": hp["value"]
            for hp in hubbard_parameters
        }
        with self.hold_trait_notifications():
            self.parameters = parameters
            self.activate = True

    def reset(self):
        with self.hold_trait_notifications():
            self.activate = self.traits()["activate"].default_value
            self.eigenvalues_label = self.traits()["eigenvalues_label"].default_value
            self.parameters = {}  # TODO default parameters
            self.eigenvalues = self._get_default_eigenvalues()

    def _get_default_parameters(self):
        return deepcopy(self._default_parameters)

    def _get_default_eigenvalues(self):
        return deepcopy(self._default_eigenvalues)

    def _get_labels(self):
        """Get a list of labels for the Hubbard widget.

        Returns:
            labels (list):
                A list of labels in the format "{kind} - {manifold}".
        """
        kind_list = self.input_structure.get_kind_names()
        hubbard_manifold_list = [
            self._get_manifold(Element(kind.symbol))
            for kind in self.input_structure.kinds
        ]
        labels = [
            f"{kind} - {manifold}"
            for kind, manifold in zip(kind_list, hubbard_manifold_list)
        ]
        return labels

    def _get_manifold(self, element):
        """Get the Hubbard manifold for a given element.

        Parameters:
            element (Element):
                The element for which to determine the Hubbard manifold.

        Returns:
            hubbard_manifold (str):
                The Hubbard manifold for the given element.
        """
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


class PseudosModel(tl.HasTraits):
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
    library = tl.Unicode(
        " ".join(
            [
                DEFAULT["advanced"]["pseudo_family"]["library"],
                DEFAULT["advanced"]["pseudo_family"]["accuracy"],
            ]
        )
    )
    functional = tl.Unicode(DEFAULT["advanced"]["pseudo_family"]["functional"])
    cutoffs = tl.List(
        trait=tl.List(tl.Float),  # [[ecutwfc values], [ecutrho values]]
        default_value=[[0.0], [0.0]],
    )
    ecutwfc = tl.Float()
    ecutrho = tl.Float()
    status_message = tl.Unicode()

    _default_dictionary = {}
    _default_cutoffs = []

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

    def set_defaults_from_structure(self):
        if self.input_structure is None:
            self._default_dictionary = {}
            self._default_cutoffs = [[0.0], [0.0]]
        else:
            self.update_pseudos()
            self.update_cutoffs()

    def update_family(self):
        library, accuracy = self.library.split()
        functional = self.functional
        # XXX (jusong.yu): a validator is needed to check the family string is consistent with the list of pseudo families defined in the setup_pseudos.py
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

        self.family = pseudo_family_string

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

        with self.hold_trait_notifications():
            self.library = f"{pseudo_family.library} {pseudo_family.accuracy}"
            self.functional = pseudo_family.functional

    def update_pseudos(self):
        try:
            pseudo_family = self._get_pseudo_family()
            pseudos = pseudo_family.get_pseudos(structure=self.input_structure)
        except ValueError as exception:
            self._status_message = f"""
                <div class='alert alert-danger'>
                    ERROR: {exception!s}
                </div>
            """
            return

        self._default_dictionary = {
            kind: pseudo.uuid for kind, pseudo in pseudos.items()
        }
        self.dictionary = self._get_default_dictionary()

    def update_cutoffs(self):
        """Update wavefunction and density cutoffs from pseudo family."""
        try:
            pseudo_family = self._get_pseudo_family()
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
            kind_names = (
                self.input_structure.get_kind_names() if self.input_structure else []
            )

        ecutwfc_list = []
        ecutrho_list = []
        for kind in kind_names:
            cutoff = cutoff_dict.get(kind, {})
            ecutrho, ecutwfc = (
                U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                for v in cutoff.values()
            )
            ecutwfc_list.append(ecutwfc)
            ecutrho_list.append(ecutrho)

        self._default_cutoffs = [ecutwfc_list or [0.0], ecutrho_list or [0.0]]
        self.cutoffs = self._get_default_cutoffs()

    def reset(self):
        with self.hold_trait_notifications():
            self.dictionary = self._get_default_dictionary()
            self.cutoffs = self._get_default_cutoffs()
            self.family = self.traits()["family"].default_value
            self.library = self.traits()["library"].default_value
            self.functional = self.traits()["functional"].default_value
            self.status_message = ""

    def _get_pseudo_family(self):
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
        return deepcopy(self._default_dictionary)

    def _get_default_cutoffs(self):
        return deepcopy(self._default_cutoffs)


class AdvancedModel(PanelModel):
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
    kpoints_distance = tl.Float(0.0)
    mesh_grid = tl.Unicode("")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.smearing = SmearingModel()
        ipw.dlink(
            (self, "protocol"),
            (self.smearing, "protocol"),
        )
        ipw.dlink(
            (self, "override"),
            (self.smearing, "override"),
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
            (self, "override"),
            (self.magnetization, "override"),
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

        self.pseudos = PseudosModel()
        ipw.dlink(
            (self, "input_structure"),
            (self.pseudos, "input_structure"),
        )
        ipw.dlink(
            (self, "protocol"),
            (self.smearing, "protocol"),
        )
        ipw.dlink(
            (self, "spin_orbit"),
            (self.pseudos, "spin_orbit"),
        )
        ipw.dlink(
            (self, "override"),
            (self.pseudos, "override"),
        )

        self.update()

    def update(self):
        parameters = PwBaseWorkChain.get_protocol_inputs(self.protocol)

        self.kpoints_distance = parameters["kpoints_distance"]

        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        etot_value = num_atoms * parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        self._set_value_and_step("etot_conv_thr", etot_value)

        scf_value = num_atoms * parameters["meta_parameters"]["conv_thr_per_atom"]
        self._set_value_and_step("scf_conv_thr", scf_value)

        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        self._set_value_and_step("forc_conv_thr", forc_value)

        self.smearing.set_defaults_from_protocol()
        self.pseudos.update_family_parameters()
        self.update_kpoints_mesh()

    def update_kpoints_mesh(self, _=None):
        if self.input_structure is None:
            self.mesh_grid = ""
        elif self.kpoints_distance > 0:
            # To avoid creating an aiida node every time we change the kpoints_distance,
            # we use the function itself instead of the decorated calcfunction.
            mesh = create_kpoints_from_distance.process_class._func(
                self.input_structure,
                orm.Float(self.kpoints_distance),
                orm.Bool(False),
            )
            self.mesh_grid = "Mesh " + str(mesh.get_kpoints_mesh()[0])
        else:
            self.mesh_grid = "Please select a number higher than 0.0"

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
                    },
                }
            },
            "clean_workdir": self.clean_workdir,
            "pseudo_family": self.pseudos.family,
            "kpoints_distance": self.kpoints_distance,
        }

        if self.hubbard.activate:
            parameters["hubbard_parameters"] = {"hubbard_u": self.hubbard.parameters}
            if self.hubbard.eigenvalues_label:
                parameters["pw"]["parameters"]["SYSTEM"].update(
                    {"starting_ns_eigenvalue": self.hubbard.eigenvalues}
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
        """Set the panel value from the given parameters."""

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
        #
        self.kpoints_distance = parameters.get("kpoints_distance", 0.15)
        if parameters.get("pw") is not None:
            system = parameters["pw"]["parameters"]["SYSTEM"]
            if "degauss" in system:
                self.smearing.degauss = system["degauss"]
            if "smearing" in system:
                self.smearing.type = system["smearing"]
            self.total_charge = parameters["pw"]["parameters"]["SYSTEM"].get(
                "tot_charge", 0
            )
            self.spin_orbit = "soc" if "lspinorb" in system else "wo_soc"
            # van der waals correction
            self.van_der_waals = self.dftd3_version.get(
                system.get("dftd3_version"),
                parameters["pw"]["parameters"]["SYSTEM"].get("vdw_corr", "none"),
            )

            # convergence threshold setting
            self.forc_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("forc_conv_thr", 0.0)
            )
            self.etot_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("CONTROL", {})
                .get("etot_conv_thr", 0.0)
            )
            self.scf_conv_thr = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("ELECTRONS", {})
                .get("conv_thr", 0.0)
            )

        # Logic to set the magnetization
        if magnetic_moments := parameters.get("initial_magnetic_moments"):
            if isinstance(magnetic_moments, list):
                magnetic_moments = {
                    kind: magnetic_moments[i]
                    for i, kind in enumerate(self.input_structure.get_kind_names())
                }
            self.magnetization.moments = magnetic_moments

        if "tot_magnetization" in parameters["pw"]["parameters"]["SYSTEM"]:
            self.magnetization.type = "tot_magnetization"

        if parameters.get("hubbard_parameters"):
            self.hubbard.activate = True
            self.hubbard.parameters = parameters["hubbard_parameters"]["hubbard_u"]
            starting_ns_eigenvalue = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("SYSTEM", {})
                .get("starting_ns_eigenvalue")
            )
            if starting_ns_eigenvalue is not None:
                self.hubbard.eigenvalues_label = True
                self.hubbard.eigenvalues = starting_ns_eigenvalue

    def reset(self):
        with self.hold_trait_notifications():
            self.override = self.traits()["override"].default_value
            self.total_charge = self.traits()["total_charge"].default_value
            self.van_der_waals = self.traits()["van_der_waals"].default_value
            self.spin_orbit = self.traits()["spin_orbit"].default_value
            self.forc_conv_thr = self.traits()["forc_conv_thr"].default_value
            self.forc_conv_thr_step = self.traits()["forc_conv_thr_step"].default_value
            self.etot_conv_thr = self.traits()["etot_conv_thr"].default_value
            self.etot_conv_thr_step = self.traits()["etot_conv_thr_step"].default_value
            self.scf_conv_thr = self.traits()["scf_conv_thr"].default_value
            self.scf_conv_thr_step = self.traits()["scf_conv_thr_step"].default_value
            self.kpoints_distance = self.traits()["kpoints_distance"].default_value
            self.update()

    def _set_value_and_step(self, attribute, value):
        """Sets the value and step size.

        Parameters:
            attribute (str):
                The attribute whose values are to be set (e.g., self.etot_conv_thr).
            value (float):
                The numerical value to set.
        """
        setattr(self, attribute, value)
        if value != 0:
            order_of_magnitude = np.floor(np.log10(abs(value)))
            setattr(self, f"{attribute}_step", 10 ** (order_of_magnitude - 1))
        else:
            setattr(self, f"{attribute}_step", 0.1)


class ConfigurationModel(tl.HasTraits):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    configuration_parameters = tl.Dict(
        key_trait=tl.Unicode(),  # parameter name
        value_trait=tl.Dict(),  # parameter configuration
        default_value={},
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.workchain = WorkChainModel()
        self.advanced = AdvancedModel()

        ipw.dlink(
            (self, "input_structure"),
            (self.advanced, "input_structure"),
        )
        ipw.dlink(
            (self.workchain, "protocol"),
            (self.advanced, "protocol"),
        )
        ipw.dlink(
            (self.workchain, "spin_type"),
            (self.advanced, "spin_type"),
        )
        ipw.dlink(
            (self.workchain, "electronic_type"),
            (self.advanced, "electronic_type"),
        )

        self._models = {
            "workchain": self.workchain,
            "advanced": self.advanced,
        }

    def add_model(self, identifier, model):
        self._models[identifier] = model

    def remove_model(self, identifier):
        if identifier in self._models:
            del self._models[identifier]

    def get_model(self, identifier):
        return self._models.get(identifier)

    def get_model_state(self):
        parameters = {
            identifier: model.get_model_state()
            for identifier, model in self._models.items()
        }
        # TODO necessary?
        parameters["workchain"].update({"properties": self._get_properties()})
        return parameters

    def set_model_state(self, parameters):
        # TODO check logic
        with self.hold_trait_notifications():
            for identifier, model in self._models.items():
                if parameters.get(identifier):
                    model.set_model_state(parameters[identifier])
            properties = parameters.get("properties", [])
            for identifier, model in self._models.items():
                if identifier in properties:
                    model.include_plugin = True
                else:
                    model.include_plugin = False

    def reset(self):
        with self.hold_trait_notifications():
            for model in self._models.values():
                model.reset()
            self.configuration_parameters = {}

    def _get_properties(self):
        properties = []
        run_bands = False
        run_pdos = False
        for identifier, model in self._models.items():
            if model.include_plugin:
                properties.append(identifier)
            if identifier == "bands":
                run_bands = True
            elif identifier == "pdos":
                run_bands = True

        if RelaxType(self.workchain.relax_type) is not RelaxType.NONE or not (
            run_bands or run_pdos
        ):
            properties.append("relax")
        return properties
