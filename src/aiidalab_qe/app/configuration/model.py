import traitlets as tl

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS as DEFAULT
from aiidalab_widgets_base import WizardAppWidgetStep


class MagnetizationSettings(tl.HasTraits):
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

    # Smearing
    degauss = tl.Float()
    smearing = tl.Unicode("gaussian")

    magnetization = MagnetizationSettings()

    # Hubbard
    activate_hubbard = tl.Bool(False)
    eigenvalues_label = tl.Bool(False)
    input_labels = tl.List(
        trait=tl.Unicode(),  # element symbol
        default_value=[],
    )
    hubbard_dict = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Float(),  # U value
        default_value={},
    )
    eigenvalues_list = tl.List(tl.List(), [])

    # Pseudopotential
    pseudos = tl.Dict(
        key_trait=tl.Unicode(),  # element symbol
        value_trait=tl.Unicode(),  # pseudopotential node uuid
        default_value={},
    )
    pseudo_family = tl.Unicode(
        "/".join(
            [
                DEFAULT["advanced"]["pseudo_family"]["library"],
                str(DEFAULT["advanced"]["pseudo_family"]["version"]),
                DEFAULT["advanced"]["pseudo_family"]["functional"],
                DEFAULT["advanced"]["pseudo_family"]["accuracy"],
            ]
        )
    )
    pseudo_library = tl.Unicode(
        " ".join(
            [
                DEFAULT["advanced"]["pseudo_family"]["library"],
                DEFAULT["advanced"]["pseudo_family"]["accuracy"],
            ]
        )
    )
    pseudo_override = tl.Bool(False)
    dft_functional = tl.Unicode(DEFAULT["advanced"]["pseudo_family"]["functional"])
    ecutwfc = tl.Float(0.0)
    ecutrho = tl.Float(0.0)


config_model = ConfigurationModel()
