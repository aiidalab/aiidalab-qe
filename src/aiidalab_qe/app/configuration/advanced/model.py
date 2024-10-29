from __future__ import annotations

import typing as t

import ipywidgets as ipw
import numpy as np
import traitlets as tl

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.panel import SettingsModel
from aiidalab_qe.setup.pseudos import PseudoFamily

if t.TYPE_CHECKING:
    from .hubbard.hubbard import HubbardModel
    from .magnetization import MagnetizationModel
    from .pseudos.pseudos import PseudosModel
    from .smearing import SmearingModel
    from .subsettings import AdvancedSubModel

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

        self._models: dict[str, AdvancedSubModel] = {}

        self._defaults = {
            "forc_conv_thr": self.traits()["forc_conv_thr"].default_value,
            "forc_conv_thr_step": self.traits()["forc_conv_thr_step"].default_value,
            "etot_conv_thr": self.traits()["etot_conv_thr"].default_value,
            "etot_conv_thr_step": self.traits()["etot_conv_thr_step"].default_value,
            "scf_conv_thr": self.traits()["scf_conv_thr"].default_value,
            "scf_conv_thr_step": self.traits()["scf_conv_thr_step"].default_value,
            "kpoints_distance": self.traits()["kpoints_distance"].default_value,
        }

    def update(self, specific=""):
        with self.hold_trait_notifications():
            self._update_defaults(specific)
            self.forc_conv_thr = self._defaults["forc_conv_thr"]
            self.forc_conv_thr_step = self._defaults["forc_conv_thr_step"]
            self.etot_conv_thr = self._defaults["etot_conv_thr"]
            self.etot_conv_thr_step = self._defaults["etot_conv_thr_step"]
            self.scf_conv_thr = self._defaults["scf_conv_thr"]
            self.scf_conv_thr_step = self._defaults["scf_conv_thr_step"]
            self.kpoints_distance = self._defaults["kpoints_distance"]

    def add_model(self, identifier, model):
        self._models[identifier] = model
        self._link_model(model)

    def get_models(self):
        return self._models.items()

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
            "kpoints_distance": self.kpoints_distance,
        }

        hubbard: HubbardModel = self._get_model("hubbard")  # type: ignore
        if hubbard.is_active:
            parameters["hubbard_parameters"] = {"hubbard_u": hubbard.parameters}
            if hubbard.has_eigenvalues:
                parameters["pw"]["parameters"]["SYSTEM"] |= {
                    "starting_ns_eigenvalue": hubbard.get_active_eigenvalues()
                }

        pseudos: PseudosModel = self._get_model("pseudos")  # type: ignore
        parameters["pseudo_family"] = pseudos.family
        if pseudos.dictionary:
            parameters["pw"]["pseudos"] = pseudos.dictionary
            parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = pseudos.ecutwfc
            parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"] = pseudos.ecutrho

        if self.van_der_waals in ["none", "ts-vdw"]:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = self.van_der_waals
        else:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = "dft-d3"
            parameters["pw"]["parameters"]["SYSTEM"]["dftd3_version"] = (
                self.dftd3_version[self.van_der_waals]
            )

        smearing: SmearingModel = self._get_model("smearing")  # type: ignore
        magnetization: MagnetizationModel = self._get_model("magnetization")  # type: ignore
        if self.spin_type == "collinear":
            parameters["initial_magnetic_moments"] = magnetization.moments
        if self.electronic_type == "metal":
            # smearing type setting
            parameters["pw"]["parameters"]["SYSTEM"]["smearing"] = smearing.type
            # smearing degauss setting
            parameters["pw"]["parameters"]["SYSTEM"]["degauss"] = smearing.degauss

        # Set tot_magnetization for collinear simulations.
        if self.spin_type == "collinear":
            # Conditions for metallic systems.
            # Select the magnetization type and set the value if override is True
            if self.electronic_type == "metal" and self.override:
                if magnetization.type == "tot_magnetization":
                    parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                        magnetization.total
                    )
                else:
                    parameters["initial_magnetic_moments"] = magnetization.moments
            # Conditions for insulator systems. Default value is 0.0
            elif self.electronic_type == "insulator":
                parameters["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                    magnetization.total
                )

        # Spin-Orbit calculation
        if self.spin_orbit == "soc":
            parameters["pw"]["parameters"]["SYSTEM"]["lspinorb"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["noncolin"] = True
            parameters["pw"]["parameters"]["SYSTEM"]["nspin"] = 4

        return parameters

    def set_model_state(self, parameters):
        pseudos: PseudosModel = self._get_model("pseudos")  # type: ignore
        if "pseudo_family" in parameters:
            pseudo_family = PseudoFamily.from_string(parameters["pseudo_family"])
            library = pseudo_family.library
            accuracy = pseudo_family.accuracy
            pseudos.library = f"{library} {accuracy}"
            pseudos.functional = pseudo_family.functional

        if "pseudos" in parameters["pw"]:
            pseudos.dictionary = parameters["pw"]["pseudos"]
            pseudos.ecutwfc = parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
            pseudos.ecutrho = parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"]

        self.kpoints_distance = parameters.get("kpoints_distance", 0.15)

        if (pw_parameters := parameters.get("pw", {}).get("parameters")) is not None:
            self._set_pw_parameters(pw_parameters)

        magnetization: MagnetizationModel = self._get_model("magnetization")  # type: ignore
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
            magnetization.moments = magnetic_moments

        hubbard: HubbardModel = self._get_model("hubbard")  # type: ignore
        if parameters.get("hubbard_parameters"):
            hubbard.is_active = True
            hubbard.parameters = parameters["hubbard_parameters"]["hubbard_u"]
            starting_ns_eigenvalue = (
                parameters.get("pw", {})
                .get("parameters", {})
                .get("SYSTEM", {})
                .get("starting_ns_eigenvalue")
            )
            if starting_ns_eigenvalue is not None:
                hubbard.has_eigenvalues = True
                hubbard.eigenvalues = starting_ns_eigenvalue

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

    def _link_model(self, model):
        ipw.dlink(
            (self, "override"),
            (model, "override"),
        )
        model.observe(
            self.unconfirm,
            tl.All,
        )
        for trait in model.dependencies:
            ipw.dlink(
                (self, trait),
                (model, trait),
            )

    def _update_defaults(self, specific=""):
        if not specific or specific != "mesh":
            parameters = PwBaseWorkChain.get_protocol_inputs(self.protocol)
            self._update_kpoints_distance(parameters)

        self._update_kpoints_mesh()

        if not specific or specific == "protocol":
            self._update_thresholds(parameters)

    def _update_kpoints_mesh(self, _=None):
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

    def _update_kpoints_distance(self, parameters):
        kpoints_distance = parameters["kpoints_distance"] if self.has_pbc else 100.0
        self._defaults["kpoints_distance"] = kpoints_distance

    def _update_thresholds(self, parameters):
        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        etot_value = num_atoms * parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        self._set_value_and_step("etot_conv_thr", etot_value)

        scf_value = num_atoms * parameters["meta_parameters"]["conv_thr_per_atom"]
        self._set_value_and_step("scf_conv_thr", scf_value)

        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        self._set_value_and_step("forc_conv_thr", forc_value)

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

        smearing: SmearingModel = self._get_model("smearing")  # type: ignore
        if "degauss" in system_params:
            smearing.degauss = system_params["degauss"]

        if "smearing" in system_params:
            smearing.type = system_params["smearing"]

        magnetization: MagnetizationModel = self._get_model("magnetization")  # type: ignore
        if "tot_magnetization" in system_params:
            magnetization.type = "tot_magnetization"

    def _get_model(self, identifier) -> AdvancedSubModel:
        if identifier in self._models:
            return self._models[identifier]
        raise ValueError(f"Model with identifier '{identifier}' not found.")
