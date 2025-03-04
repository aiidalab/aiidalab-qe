from __future__ import annotations

import typing as t

import ipywidgets as ipw
import numpy as np
import traitlets as tl

from aiida import orm
from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import HasInputStructure, HasModels
from aiidalab_qe.common.panel import ConfigurationSettingsModel
from aiidalab_qe.setup.pseudos import PseudoFamily

from .subsettings import AdvancedCalculationSubSettingsModel

if t.TYPE_CHECKING:
    from .hubbard.hubbard import HubbardConfigurationSettingsModel
    from .magnetization import MagnetizationConfigurationSettingsModel
    from .pseudos.pseudos import PseudosConfigurationSettingsModel
    from .smearing import SmearingConfigurationSettingsModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class AdvancedConfigurationSettingsModel(
    ConfigurationSettingsModel,
    HasModels[AdvancedCalculationSubSettingsModel],
    HasInputStructure,
):
    title = "Advanced settings"
    identifier = "advanced"

    dependencies = [
        "input_structure",
        "workchain.protocol",
        "workchain.spin_type",
        "workchain.electronic_type",
        "workchain.spin_orbit",
    ]

    protocol = tl.Unicode()
    spin_type = tl.Unicode()
    electronic_type = tl.Unicode()
    spin_orbit = tl.Unicode()

    clean_workdir = tl.Bool(DEFAULT["advanced"]["clean_workdir"])
    total_charge = tl.Float(DEFAULT["advanced"]["tot_charge"])
    van_der_waals_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["None", "none"],
            ["Grimme-D3", "dft-d3"],
            ["Grimme-D3BJ", "dft-d3bj"],
            ["Grimme-D3M", "dft-d3m"],
            ["Grimme-D3MBJ", "dft-d3mbj"],
            ["Tkatchenko-Scheffler", "ts-vdw"],
        ],
    )
    van_der_waals = tl.Unicode(DEFAULT["advanced"]["vdw_corr"])
    forc_conv_thr = tl.Float(0.0)
    forc_conv_thr_step = tl.Float(1e-4)
    etot_conv_thr = tl.Float(0.0)
    etot_conv_thr_step = tl.Float(1e-5)
    scf_conv_thr = tl.Float(0.0)
    scf_conv_thr_step = tl.Float(1e-10)
    electron_maxstep = tl.Int(80)
    optimization_maxsteps = tl.Int(50)

    kpoints_distance = tl.Float(0.0)
    mesh_grid = tl.Unicode("")

    include = True

    dftd3_version = {
        "dft-d3": 3,
        "dft-d3bj": 4,
        "dft-d3m": 5,
        "dft-d3mbj": 6,
    }

    def update(self, specific=""):
        with self.hold_trait_notifications():
            if not specific or specific != "mesh":
                parameters = PwBaseWorkChain.get_protocol_inputs(self.protocol)
                self._update_kpoints_distance(parameters)
                if specific == "protocol":
                    self._update_thresholds(parameters)
            self._update_kpoints_mesh()

    def get_model_state(self):
        num_atoms = len(self.input_structure.sites) if self.input_structure else 1
        parameters = {
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "tot_charge": self.total_charge,
                    },
                    "CONTROL": {
                        "forc_conv_thr": self.forc_conv_thr,
                        "etot_conv_thr": self.etot_conv_thr * num_atoms,
                    },
                    "ELECTRONS": {
                        "conv_thr": self.scf_conv_thr * num_atoms,
                        "electron_maxstep": self.electron_maxstep,
                    },
                }
            },
            "clean_workdir": self.clean_workdir,
            "kpoints_distance": self.kpoints_distance,
            "optimization_maxsteps": self.optimization_maxsteps,
        }

        hubbard: HubbardConfigurationSettingsModel = self.get_model("hubbard")  # type: ignore
        if hubbard.is_active:
            parameters["hubbard_parameters"] = {
                "hubbard_u": {
                    label: value for label, value in hubbard.parameters.items() if value
                }
            }
            if hubbard.has_eigenvalues:
                parameters["pw"]["parameters"]["SYSTEM"] |= {
                    "starting_ns_eigenvalue": hubbard.get_active_eigenvalues()
                }

        pseudos: PseudosConfigurationSettingsModel = self.get_model("pseudos")  # type: ignore
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

        smearing: SmearingConfigurationSettingsModel = self.get_model("smearing")  # type: ignore
        if self.electronic_type == "metal":
            # smearing type setting
            parameters["pw"]["parameters"]["SYSTEM"]["smearing"] = smearing.type
            # smearing degauss setting
            parameters["pw"]["parameters"]["SYSTEM"]["degauss"] = smearing.degauss

        magnetization: MagnetizationConfigurationSettingsModel = self.get_model(
            "magnetization"
        )  # type: ignore
        if self.spin_type == "collinear":
            parameters["initial_magnetic_moments"] = magnetization.moments

        # Set tot_magnetization for collinear simulations.
        if self.spin_type == "collinear":
            # Conditions for metallic systems.
            # Select the magnetization type and set the value
            if self.electronic_type == "metal":
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
        pseudos: PseudosConfigurationSettingsModel = self.get_model("pseudos")  # type: ignore
        if "pseudo_family" in parameters:
            pseudo_family = PseudoFamily.from_string(parameters["pseudo_family"])
            library = pseudo_family.library
            accuracy = pseudo_family.accuracy
            pseudos.library = f"{library} {accuracy}"
            pseudos.functional = pseudo_family.functional
            pseudos.family = parameters["pseudo_family"]

        if "pseudos" in parameters["pw"]:
            pseudos.dictionary = parameters["pw"]["pseudos"]
            pseudos.ecutwfc = parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
            pseudos.ecutrho = parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"]

        self.kpoints_distance = parameters.get("kpoints_distance", 0.15)
        self.optimization_maxsteps = parameters.get("optimization_maxsteps", 50)

        if (pw_parameters := parameters.get("pw", {}).get("parameters")) is not None:
            self._set_pw_parameters(pw_parameters)

        magnetization: MagnetizationConfigurationSettingsModel = self.get_model(
            "magnetization"
        )  # type: ignore
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

        hubbard: HubbardConfigurationSettingsModel = self.get_model("hubbard")  # type: ignore
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
                hubbard.set_active_eigenvalues(starting_ns_eigenvalue)

    def reset(self):
        with self.hold_trait_notifications():
            self.total_charge = self._get_default("total_charge")
            self.van_der_waals = self._get_default("van_der_waals")
            self.forc_conv_thr = self._get_default("forc_conv_thr")
            self.forc_conv_thr_step = self._get_default("forc_conv_thr_step")
            self.etot_conv_thr = self._get_default("etot_conv_thr")
            self.etot_conv_thr_step = self._get_default("etot_conv_thr_step")
            self.scf_conv_thr = self._get_default("scf_conv_thr")
            self.scf_conv_thr_step = self._get_default("scf_conv_thr_step")
            self.electron_maxstep = self._get_default("electron_maxstep")
            self.kpoints_distance = self._get_default("kpoints_distance")
            self.optimization_maxsteps = self._get_default("optimization_maxsteps")

    def _get_default(self, trait):
        return self._defaults.get(trait, self.traits()[trait].default_value)

    def _link_model(self, model: AdvancedCalculationSubSettingsModel):
        ipw.dlink(
            (self, "loaded_from_process"),
            (model, "loaded_from_process"),
        )
        model.observe(
            self._on_any_change,
            tl.All,
        )
        super()._link_model(model)

    def _update_kpoints_mesh(self, _=None):
        if not self.has_structure:
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
        self.kpoints_distance = self._defaults["kpoints_distance"]

    def _update_thresholds(self, parameters):
        etot_value = parameters["meta_parameters"]["etot_conv_thr_per_atom"]
        self._set_value_and_step("etot_conv_thr", etot_value)
        self.etot_conv_thr = self._defaults["etot_conv_thr"]
        self.etot_conv_thr_step = self._defaults["etot_conv_thr_step"]

        scf_value = parameters["meta_parameters"]["conv_thr_per_atom"]
        self._set_value_and_step("scf_conv_thr", scf_value)
        self.scf_conv_thr = self._defaults["scf_conv_thr"]
        self.scf_conv_thr_step = self._defaults["scf_conv_thr_step"]

        forc_value = parameters["pw"]["parameters"]["CONTROL"]["forc_conv_thr"]
        self._set_value_and_step("forc_conv_thr", forc_value)
        self.forc_conv_thr = self._defaults["forc_conv_thr"]
        self.forc_conv_thr_step = self._defaults["forc_conv_thr_step"]

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

        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        self.forc_conv_thr = control_params.get("forc_conv_thr", 0.0)
        self.etot_conv_thr = control_params.get("etot_conv_thr", 0.0) / num_atoms
        self.scf_conv_thr = electron_params.get("conv_thr", 0.0) / num_atoms
        self.electron_maxstep = electron_params.get("electron_maxstep", 80)

        self.total_charge = system_params.get("tot_charge", 0)
        self.spin_orbit = "soc" if "lspinorb" in system_params else "wo_soc"

        self.van_der_waals = self.dftd3_version.get(
            system_params.get("dftd3_version"),
            system_params.get("vdw_corr", "none"),
        )

        smearing: SmearingConfigurationSettingsModel = self.get_model("smearing")  # type: ignore
        if "degauss" in system_params:
            smearing.degauss = system_params["degauss"]

        if "smearing" in system_params:
            smearing.type = system_params["smearing"]

        magnetization: MagnetizationConfigurationSettingsModel = self.get_model(
            "magnetization"
        )  # type: ignore
        if "tot_magnetization" in system_params:
            magnetization.type = "tot_magnetization"
            magnetization.total = system_params["tot_magnetization"]
