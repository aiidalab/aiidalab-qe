from __future__ import annotations

import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.mixins import HasInputStructure, HasModels
from aiidalab_qe.common.panel import ConfigurationSettingsModel
from aiidalab_qe.setup.pseudos import PseudoFamily
from aiidalab_qe.utils import get_pseudo_info

from .convergence import ConvergenceConfigurationSettingsModel
from .general import GeneralConfigurationSettingsModel
from .hubbard import HubbardConfigurationSettingsModel
from .magnetization import MagnetizationConfigurationSettingsModel
from .pseudos import PseudosConfigurationSettingsModel
from .smearing import SmearingConfigurationSettingsModel
from .subsettings import AdvancedCalculationSubSettingsModel


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

    include = True

    def get_model_state(self):
        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        general = t.cast(
            GeneralConfigurationSettingsModel,
            self.get_model("general"),
        )
        convergence = t.cast(
            ConvergenceConfigurationSettingsModel,
            self.get_model("convergence"),
        )
        parameters = {
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "tot_charge": general.total_charge,
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
            "clean_workdir": general.clean_workdir,
            "kpoints_distance": convergence.kpoints_distance,
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

        if general.van_der_waals in ["none", "ts-vdw"]:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = general.van_der_waals
        else:
            parameters["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = "dft-d3"
            parameters["pw"]["parameters"]["SYSTEM"]["dftd3_version"] = (
                general.dftd3_version[general.van_der_waals]
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
            library = f"{pseudo_family.library} {pseudo_family.accuracy}"
            if relativistic := pseudo_family.relativistic:
                library += f" ({relativistic})"
            pseudos.functional = pseudo_family.functional
            pseudos.library = library
            pseudos.family = pseudo_family_string
        else:
            pp_uuid = next(iter(parameters["pw"]["pseudos"].values()))
            pseudo_info = get_pseudo_info(pp_uuid)
            pseudos.functional = pseudo_info["functional"]
            pseudos.library = None
            pseudos.family = None
            pseudos.show_upload_warning = True

        pseudos.functionals = [pseudos.functional] * len(pseudos.dictionary)

        if "pseudos" in parameters["pw"]:
            pseudos.dictionary = parameters["pw"]["pseudos"]
            pseudos.ecutwfc = parameters["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
            pseudos.ecutrho = parameters["pw"]["parameters"]["SYSTEM"]["ecutrho"]

        convergence = t.cast(
            ConvergenceConfigurationSettingsModel,
            self.get_model("convergence"),
        )
        convergence.optimization_maxsteps = parameters.get("optimization_maxsteps", 50)
        convergence.kpoints_distance = parameters.get("kpoints_distance", 0.15)

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

    def _set_pw_parameters(self, pw_parameters):
        system_params: dict = pw_parameters.get("SYSTEM", {})
        control_params: dict = pw_parameters.get("CONTROL", {})
        electron_params: dict = pw_parameters.get("ELECTRONS", {})

        num_atoms = len(self.input_structure.sites) if self.input_structure else 1

        general = t.cast(
            GeneralConfigurationSettingsModel,
            self.get_model("general"),
        )
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

        general.total_charge = system_params.get("tot_charge", 0)

        # NOTE: this is here for backwards compatability, as SYSTEM parameters
        # are part of the advanced input namespace.
        self.spin_orbit = "soc" if "lspinorb" in system_params else "wo_soc"

        general.van_der_waals = general.dftd3_version.get(
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
