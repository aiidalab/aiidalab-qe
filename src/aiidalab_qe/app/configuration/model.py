import traitlets as tl
from pymatgen.core.periodic_table import Element

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS as DEFAULT
from aiidalab_widgets_base import WizardAppWidgetStep


class SmearingModel(tl.HasTraits):
    type = tl.Unicode()
    degauss = tl.Float()

    def set_defaults(self, default_protocol):
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

    def reset(self):
        self.type = self.traits()["type"].default_value
        self.total = self.traits()["total"].default_value
        self.moments = {}


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

    def set_defaults(self, change):
        if (input_structure := change["new"]) is None:
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
        self.eigenvalues = self._default_eigenvalues
        self.needs_eigenvalues_widget = len(self.elements) > 0

    def reset(self):
        self.activate = self.traits()["activate"].default_value
        self.eigenvalues_label = self.traits()["eigenvalues_label"].default_value
        self.parameters = {}
        self.eigenvalues = self._default_eigenvalues


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
    total_charge = tl.Int(DEFAULT["advanced"]["tot_charge"])
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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.smearing.set_defaults(self.traits()["protocol"].default_value)
        self.observe(self.hubbard.set_defaults, "input_structure")


config_model = ConfigurationModel()
