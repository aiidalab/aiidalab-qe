from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl
from aiida_pseudo.common.units import U
from pymatgen.core.periodic_table import Element

from aiida import orm
from aiida.common import exceptions
from aiida.plugins import GroupFactory
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_widgets_base import WizardAppWidgetStep

SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SmearingModel(tl.HasTraits):
    type = tl.Unicode()
    degauss = tl.Float()

    def set_defaults_for_protocol(self, default_protocol):
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

    @tl.default("type")
    def _default_type(self):
        return self._defaults["type"]

    @tl.default("degauss")
    def _default_degauss(self):
        return self._defaults["degauss"]

    def reset(self):
        self.type = self._default_type()
        self.degauss = self._default_degauss()


class MagnetizationModel(tl.HasTraits):
    type = tl.Unicode("starting_magnetization")
    total = tl.Float(0.0)
    moments = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # magnetic moment
        default_value={},
    )

    def set_defaults_for_structure(self, input_structure):
        if input_structure is None:
            self._default_moments = {}
        else:
            self._default_moments = {kind.symbol: 0.0 for kind in input_structure.kinds}
        self.moments = self._get_moments()

    def _get_moments(self):
        return deepcopy(self._default_moments)

    def reset(self):
        self.type = self.traits()["type"].default_value
        self.total = self.traits()["total"].default_value
        self.moments = self._get_moments()


class HubbardModel(tl.HasTraits):
    activate = tl.Bool(False)
    eigenvalues_label = tl.Bool(False)
    parameters = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # U value
        default_value={},
    )
    eigenvalues = tl.List(
        trait=tl.List(),  # [[[[state, spin, kind, eigenvalue] # state] # spin] # kind]
        default_value=[],
    )

    def set_defaults_for_structure(self, input_structure):
        if input_structure is None:
            self.elements = []
            self._default_eigenvalues = []
        else:
            self.elements = [
                *filter(
                    lambda element: (
                        element.is_transition_metal
                        or element.is_lanthanoid
                        or element.is_actinoid
                    ),
                    [Element(kind.symbol) for kind in input_structure.kinds],
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
        self.eigenvalues = self._get_default_eigenvalues()
        self.needs_eigenvalues_widget = len(self.elements) > 0

    def _get_default_eigenvalues(self):
        return deepcopy(self._default_eigenvalues)

    def reset(self):
        self.activate = self.traits()["activate"].default_value
        self.eigenvalues_label = self.traits()["eigenvalues_label"].default_value
        self.parameters = {}  # TODO default parameters
        self.eigenvalues = self._get_default_eigenvalues()


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

    def set_defaults_for_structure(self, input_structure):
        if input_structure is None:
            self._default_dictionary = {}
            self._default_cutoffs = [[0.0], [0.0]]
        else:
            self.update_pseudos(input_structure)
            self.update_cutoffs(input_structure)

    def update_pseudos(self, input_structure):
        try:
            pseudo_family = self._get_pseudo_family()
            pseudos = pseudo_family.get_pseudos(structure=input_structure)
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

    def update_cutoffs(self, input_structure):
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
            kind_names = input_structure.get_kind_names()

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

        self._default_cutoffs = [ecutwfc_list, ecutrho_list]
        self.cutoffs = self._get_default_cutoffs()

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

    def reset(self):
        self.dictionary = self._get_default_dictionary()
        self.cutoffs = self._get_default_cutoffs()
        self.family = self.traits()["family"].default_value
        self.library = self.traits()["library"].default_value
        self.functional = self.traits()["functional"].default_value
        self.status_message = ""


class ConfigurationModel(tl.HasTraits):
    state = tl.UseEnum(WizardAppWidgetStep.State)
    prev_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    # Basic settings
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
    protocol = tl.Unicode(DEFAULT["workchain"]["protocol"])
    relax_type = tl.Unicode("positions_cell")
    spin_type = tl.Unicode(DEFAULT["workchain"]["spin_type"])
    electronic_type = tl.Unicode(DEFAULT["workchain"]["electronic_type"])

    # Advanced setting defaults
    clean_workdir = tl.Bool(False)
    override = tl.Bool(False)
    total_charge = tl.Float(DEFAULT["advanced"]["tot_charge"])
    van_der_waals = tl.Unicode(DEFAULT["advanced"]["vdw_corr"])
    spin_orbit = tl.Unicode("wo_soc")

    # Convergence
    forc_conv_thr = tl.Float(0.0)
    forc_conv_thr_step = tl.Float(1e-4)
    etot_conv_thr = tl.Float(0.0)
    etot_conv_thr_step = tl.Float(1e-5)
    scf_conv_thr = tl.Float(0.0)
    scf_conv_thr_step = tl.Float(1e-10)

    # Kpoints
    kpoints_distance = tl.Float(0.0)
    mesh_grid = tl.Unicode()

    smearing = SmearingModel()
    magnetization = MagnetizationModel()
    hubbard = HubbardModel()
    pseudos = PseudosModel()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.smearing.set_defaults_for_protocol(self.traits()["protocol"].default_value)


config_model = ConfigurationModel()
