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
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION, PseudoFamily
from aiidalab_widgets_base import WizardAppWidgetStep

SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class BasicModel(tl.HasTraits):
    protocol = tl.Unicode(DEFAULT["workchain"]["protocol"])
    relax_type = tl.Unicode("positions_cell")
    spin_type = tl.Unicode(DEFAULT["workchain"]["spin_type"])
    electronic_type = tl.Unicode(DEFAULT["workchain"]["electronic_type"])

    def reset(self):
        with self.hold_trait_notifications():
            self.protocol = self.traits()["protocol"].default_value
            self.relax_type = self.traits()["relax_type"].default_value
            self.spin_type = self.traits()["spin_type"].default_value
            self.electronic_type = self.traits()["electronic_type"].default_value


class AdvancedModel(tl.HasTraits):
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

    def update_kpoints_mesh(self, structure, distance=None):
        if structure is None:
            self.mesh_grid = ""
        else:
            distance = self.kpoints_distance if distance is None else distance
            if distance > 0:
                # To avoid creating an aiida node every time we change the kpoints_distance,
                # we use the function itself instead of the decorated calcfunction.
                mesh = create_kpoints_from_distance.process_class._func(
                    structure,
                    orm.Float(distance),
                    orm.Bool(False),
                )
                self.mesh_grid = "Mesh " + str(mesh.get_kpoints_mesh()[0])
            else:
                self.mesh_grid = "Please select a number higher than 0.0"

    def set_value_and_step(self, attribute, value):
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

    def reset(self):
        with self.hold_trait_notifications():
            self.clean_workdir = self.traits()["clean_workdir"].default_value
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
            self.mesh_grid = self.traits()["mesh_grid"].default_value


class SmearingModel(tl.HasTraits):
    type = tl.Unicode()
    degauss = tl.Float()

    _defaults = {}

    def set_defaults_from_protocol(self, default_protocol):
        parameters = (
            PwBaseWorkChain.get_protocol_inputs(default_protocol)
            .get("pw", {})
            .get("parameters", {})
            .get("SYSTEM", {})
        )
        self._defaults = {
            "type": parameters["smearing"],
            "degauss": parameters["degauss"],
        }

    def update_from_protocol(self, protocol):
        parameters = (
            PwBaseWorkChain.get_protocol_inputs(protocol)
            .get("pw", {})
            .get("parameters", {})
            .get("SYSTEM", {})
        )
        with self.hold_trait_notifications():
            self.type = parameters["smearing"]
            self.degauss = parameters["degauss"]

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
    type = tl.Unicode("starting_magnetization")
    total = tl.Float(0.0)
    moments = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # magnetic moment
        default_value={},
    )

    _default_moments = {}

    def set_defaults_from_structure(self, structure):
        if structure is None:
            self._default_moments = {}
        else:
            self._default_moments = {kind.symbol: 0.0 for kind in structure.kinds}
        self.moments = self._get_moments()

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self.traits()["type"].default_value
            self.total = self.traits()["total"].default_value
            self.moments = self._get_moments()

    def _get_moments(self):
        return deepcopy(self._default_moments)


class HubbardModel(tl.HasTraits):
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

    def set_defaults_from_structure(self, structure):
        if structure is None:
            self.elements = []  # TODO: rename for clarity
            self.input_labels = []  # TODO: rename for clarity
            self._default_parameters = {}
            self._default_eigenvalues = []
        else:
            self.input_labels = self._get_labels(structure)
            self._default_parameters = {label: 0.0 for label in self.input_labels}
            self.elements = [
                *filter(
                    lambda element: (
                        element.is_transition_metal
                        or element.is_lanthanoid
                        or element.is_actinoid
                    ),
                    [Element(kind.symbol) for kind in structure.kinds],
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

    def set_parameters_from_hubbard_structure(self, structure):
        hubbard_parameters = structure.hubbard.dict()["parameters"]
        sites = structure.sites
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

    def _get_labels(self, structure):
        """Get a list of labels for the Hubbard widget.

        Returns:
            labels (list):
                A list of labels in the format "{kind} - {manifold}".
        """
        kind_list = structure.get_kind_names()
        hubbard_manifold_list = [
            self._get_manifold(Element(kind.symbol)) for kind in structure.kinds
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

    def set_defaults_from_structure(self, structure):
        if structure is None:
            self._default_dictionary = {}
            self._default_cutoffs = [[0.0], [0.0]]
        else:
            self.update_pseudos(structure)
            self.update_cutoffs(structure)

    def update_family(self, spin_orbit):
        library, accuracy = self.library.split()
        functional = self.functional
        # XXX (jusong.yu): a validator is needed to check the family string is consistent with the list of pseudo families defined in the setup_pseudos.py
        if library == "PseudoDojo":
            if spin_orbit == "soc":
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

    def update_family_parameters(self, protocol, spin_orbit):
        if spin_orbit == "soc":
            if protocol in ["fast", "moderate"]:
                pseudo_family_string = "PseudoDojo/0.4/PBE/FR/standard/upf"
            else:
                pseudo_family_string = "PseudoDojo/0.4/PBE/FR/stringent/upf"
        else:
            pseudo_family_string = PwBaseWorkChain.get_protocol_inputs(protocol)[
                "pseudo_family"
            ]

        pseudo_family = PseudoFamily.from_string(pseudo_family_string)

        with self.hold_trait_notifications():
            self.library = f"{pseudo_family.library} {pseudo_family.accuracy}"
            self.functional = pseudo_family.functional

    def update_pseudos(self, structure):
        try:
            pseudo_family = self._get_pseudo_family()
            pseudos = pseudo_family.get_pseudos(structure=structure)
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

    def update_cutoffs(self, structure):
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
            kind_names = structure.get_kind_names() if structure else []

        ecutwfc_list = [0.0]
        ecutrho_list = [0.0]
        for kind in kind_names:
            cutoff = cutoff_dict.get(kind, {})
            ecutrho, ecutwfc = (
                U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                for v in cutoff.values()
            )
            ecutwfc_list.append(ecutwfc)
            ecutrho_list.append(ecutrho)

        self._default_cutoffs = [ecutwfc_list, ecutrho_list]
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


class ConfigurationModel(tl.HasTraits):
    state = tl.UseEnum(
        enum_class=WizardAppWidgetStep.State,
        default_value=WizardAppWidgetStep.State.INIT,
    )
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)
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

    basic = BasicModel()
    advanced = AdvancedModel()
    smearing = SmearingModel()
    magnetization = MagnetizationModel()
    hubbard = HubbardModel()
    pseudos = PseudosModel()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.observe(self._update_from_structure, "input_structure")
        self.update_from_protocol()

    def update_from_protocol(self, protocol=None):
        protocol = self.basic.protocol if protocol is None else protocol

        parameters = PwBaseWorkChain.get_protocol_inputs(protocol)

        self.advanced.kpoints_distance = parameters["kpoints_distance"]

        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        etot_value = num_atoms * parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        self.advanced.set_value_and_step("etot_conv_thr", etot_value)

        scf_value = num_atoms * parameters["meta_parameters"]["conv_thr_per_atom"]
        self.advanced.set_value_and_step("scf_conv_thr", scf_value)

        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        self.advanced.set_value_and_step("forc_conv_thr", forc_value)

        self.smearing.set_defaults_from_protocol(protocol)
        self.pseudos.update_family_parameters(protocol, self.advanced.spin_orbit)

    def _update_from_structure(self, change):
        self.advanced.update_kpoints_mesh(change["new"])


config_model = ConfigurationModel()
