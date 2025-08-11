from __future__ import annotations

import typing as t

import ipywidgets as ipw
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
from aiidalab_qe.utils import get_pseudo_info

from .convergence import ConvergenceConfigurationSettingsModel
from .hubbard import HubbardConfigurationSettingsModel
from .magnetization import MagnetizationConfigurationSettingsModel
from .pseudos import PseudosConfigurationSettingsModel
from .smearing import SmearingConfigurationSettingsModel
from .subsettings import AdvancedCalculationSubSettingsModel

DEFAULT = t.cast(dict, DEFAULT_PARAMETERS)


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
            self._update_kpoints_mesh()

    def get_model_state(self):
        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        convergence = t.cast(
            ConvergenceConfigurationSettingsModel,
            self.get_model("convergence"),
        )
        parameters = {
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "tot_charge": self.total_charge,
                    },
                    "CONTROL": {
                        "forc_conv_thr": convergence.forc_conv_thr,
                        "etot_conv_thr": convergence.etot_conv_thr * num_atoms,
                    },
                    "ELECTRONS": {
                        "conv_thr": convergence.scf_conv_thr * num_atoms,
                        "electron_maxstep": convergence.electron_maxstep,
                        "mixing_beta": convergence.mixing_beta,
                    },
                }
            },
            "clean_workdir": self.clean_workdir,
            "kpoints_distance": self.kpoints_distance,
            "optimization_maxsteps": convergence.optimization_maxsteps,
        }

        # Only modify if mixing mode is different than default
        if convergence.mixing_mode != "plain":
            parameters["pw"]["parameters"]["ELECTRONS"]["mixing_mode"] = (
                convergence.mixing_mode
            )

        hubbard = t.cast(
            HubbardConfigurationSettingsModel,
            self.get_model("hubbard"),
        )
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

        pseudos = t.cast(
            PseudosConfigurationSettingsModel,
            self.get_model("pseudos"),
        )
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

        smearing = t.cast(
            SmearingConfigurationSettingsModel,
            self.get_model("smearing"),
        )
        if self.electronic_type == "metal":
            # smearing type setting
            parameters["pw"]["parameters"]["SYSTEM"]["smearing"] = smearing.type
            # smearing degauss setting
            parameters["pw"]["parameters"]["SYSTEM"]["degauss"] = smearing.degauss

        magnetization = t.cast(
            MagnetizationConfigurationSettingsModel,
            self.get_model("magnetization"),
        )
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
        pseudos = t.cast(
            PseudosConfigurationSettingsModel,
            self.get_model("pseudos"),
        )
        if pseudo_family_string := parameters.get("pseudo_family"):
            pseudo_family = PseudoFamily.from_string(pseudo_family_string)
            library = pseudo_family.library
            accuracy = pseudo_family.accuracy
            pseudos.library = f"{library} {accuracy}"
            pseudos.family = pseudo_family_string
        else:
            pseudos.library = None
            pp_uuid = next(iter(parameters["pw"]["pseudos"].values()))
            pseudo_info = get_pseudo_info(pp_uuid)
            pseudos.functional = pseudo_info["functional"]
            pseudos.family = None
            pseudos.show_upload_warning = True

        if "pseudos" in parameters["pw"]:
            pseudos.dictionary = parameters["pw"]["pseudos"]
            pseudos.ecutwfc = parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
            pseudos.ecutrho = parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"]
            pseudos.functionals = [pseudos.functional] * len(pseudos.dictionary)

        convergence = t.cast(
            ConvergenceConfigurationSettingsModel,
            self.get_model("convergence"),
        )
        convergence.optimization_maxsteps = parameters.get("optimization_maxsteps", 50)
        self.kpoints_distance = parameters.get("kpoints_distance", 0.15)

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
            self.kpoints_distance = self._get_default("kpoints_distance")

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

    def _set_pw_parameters(self, pw_parameters):
        system_params: dict = pw_parameters.get("SYSTEM", {})
        control_params: dict = pw_parameters.get("CONTROL", {})
        electron_params: dict = pw_parameters.get("ELECTRONS", {})

        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        convergence = t.cast(
            ConvergenceConfigurationSettingsModel,
            self.get_model("convergence"),
        )
        convergence.forc_conv_thr = control_params.get("forc_conv_thr", 0.0)
        convergence.etot_conv_thr = control_params.get("etot_conv_thr", 0.0) / num_atoms
        convergence.scf_conv_thr = electron_params.get("conv_thr", 0.0) / num_atoms
        convergence.electron_maxstep = electron_params.get("electron_maxstep", 80)
        convergence.mixing_mode = electron_params.get("mixing_mode", "plain")
        convergence.mixing_beta = electron_params.get("mixing_beta", 0.4)

        self.total_charge = system_params.get("tot_charge", 0)
        self.spin_orbit = "soc" if "lspinorb" in system_params else "wo_soc"

        self.van_der_waals = self.dftd3_version.get(
            system_params.get("dftd3_version", ""),
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
