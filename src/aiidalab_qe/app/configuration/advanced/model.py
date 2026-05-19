from __future__ import annotations

import typing as t

import traitlets as tl

from aiidalab_qe.common.mixins import HasInputStructure, HasModels
from aiidalab_qe.common.panel import PanelModel
from aiidalab_qe.setup.pseudos import PseudoFamily
from aiidalab_qe.utils import get_pseudo_info

from .convergence import ConvergenceConfigurationSettingsModel
from .general import GeneralConfigurationSettingsModel
from .hubbard import HubbardConfigurationSettingsModel
from .magnetization import MagnetizationConfigurationSettingsModel
from .pseudos import PseudosConfigurationSettingsModel
from .smearing import SmearingConfigurationSettingsModel


class AdvancedConfigurationSettingsModel(
    PanelModel,
    HasModels[PanelModel],
    HasInputStructure,
):
    title = "Advanced settings"
    identifier = "advanced"

    dependencies = [
        "structure_uuid",
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

    def get_model_state(self) -> dict:
        num_atoms = len(self.input_structure.sites) if self.has_structure else 1

        general = t.cast(
            GeneralConfigurationSettingsModel,
            self.get_model("general"),
        )
        convergence = t.cast(
            ConvergenceConfigurationSettingsModel,
            self.get_model("convergence"),
        )
        state = {
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
            state["pw"]["parameters"]["ELECTRONS"]["mixing_mode"] = (
                convergence.mixing_mode
            )

        hubbard = t.cast(
            HubbardConfigurationSettingsModel,
            self.get_model("hubbard"),
        )
        if hubbard.is_active:
            state["hubbard_parameters"] = {
                "hubbard_u": {
                    label: value for label, value in hubbard.parameters.items() if value
                }
            }
            if hubbard.has_eigenvalues:
                state["pw"]["parameters"]["SYSTEM"] |= {
                    "starting_ns_eigenvalue": hubbard.get_active_eigenvalues()
                }

        pseudos = t.cast(
            PseudosConfigurationSettingsModel,
            self.get_model("pseudos"),
        )
        state["pseudo_family"] = pseudos.family
        if pseudos.dictionary:
            state["pw"]["pseudos"] = pseudos.dictionary
            state["pw"]["parameters"]["SYSTEM"]["ecutwfc"] = pseudos.ecutwfc
            state["pw"]["parameters"]["SYSTEM"]["ecutrho"] = pseudos.ecutrho

        if general.van_der_waals in ["none", "ts-vdw"]:
            state["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = general.van_der_waals
        else:
            state["pw"]["parameters"]["SYSTEM"]["vdw_corr"] = "dft-d3"
            state["pw"]["parameters"]["SYSTEM"]["dftd3_version"] = (
                general.dftd3_version[general.van_der_waals]
            )

        smearing = t.cast(
            SmearingConfigurationSettingsModel,
            self.get_model("smearing"),
        )
        if self.electronic_type == "metal":
            # smearing type setting
            state["pw"]["parameters"]["SYSTEM"]["smearing"] = smearing.type
            # smearing degauss setting
            state["pw"]["parameters"]["SYSTEM"]["degauss"] = smearing.degauss

        magnetization = t.cast(
            MagnetizationConfigurationSettingsModel,
            self.get_model("magnetization"),
        )
        if self.spin_type == "collinear":
            state["initial_magnetic_moments"] = magnetization.moments

        # Set tot_magnetization for collinear simulations.
        if self.spin_type == "collinear":
            # Conditions for metallic systems.
            # Select the magnetization type and set the value
            if self.electronic_type == "metal":
                if magnetization.type == "tot_magnetization":
                    state["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                        magnetization.total
                    )
                else:
                    state["initial_magnetic_moments"] = magnetization.moments
            # Conditions for insulator systems. Default value is 0.0
            elif self.electronic_type == "insulator":
                state["pw"]["parameters"]["SYSTEM"]["tot_magnetization"] = (
                    magnetization.total
                )

        # Spin-Orbit calculation
        if self.spin_orbit == "soc":
            state["pw"]["parameters"]["SYSTEM"]["lspinorb"] = True
            state["pw"]["parameters"]["SYSTEM"]["noncolin"] = True
            state["pw"]["parameters"]["SYSTEM"]["nspin"] = 4

        return state

    def set_model_state(self, state: dict):
        PW: dict = state.get("pw", {})
        PW_PARAMS: dict = PW.get("parameters", {})
        SYSTEM: dict = PW_PARAMS.get("SYSTEM", {})
        CONTROL: dict = PW_PARAMS.get("CONTROL", {})
        ELECTRONS: dict = PW_PARAMS.get("ELECTRONS", {})

        num_atoms = len(self.input_structure.sites) if self.has_structure else 1

        self.spin_orbit = "soc" if "lspinorb" in SYSTEM else "wo_soc"

        general = t.cast(
            GeneralConfigurationSettingsModel,
            self.get_model("general"),
        )
        with general.hold_trait_notifications():
            general.clean_workdir = state.get("clean_workdir", True)
            general.total_charge = SYSTEM.get("tot_charge", 0)
            general.van_der_waals = general.dftd3_version.get(
                SYSTEM.get("dftd3_version", ""),
                SYSTEM.get("vdw_corr", "none"),
            )

        convergence = t.cast(
            ConvergenceConfigurationSettingsModel,
            self.get_model("convergence"),
        )
        with convergence.hold_trait_notifications():
            convergence.optimization_maxsteps = state.get("optimization_maxsteps", 50)
            convergence.kpoints_distance = state.get("kpoints_distance", 0.15)
            convergence.forc_conv_thr = CONTROL.get("forc_conv_thr", 0.0)
            convergence.etot_conv_thr = CONTROL.get("etot_conv_thr", 0.0) / num_atoms
            convergence.scf_conv_thr = ELECTRONS.get("conv_thr", 0.0) / num_atoms
            convergence.electron_maxstep = ELECTRONS.get("electron_maxstep", 80)
            convergence.mixing_mode = ELECTRONS.get("mixing_mode", "plain")
            convergence.mixing_beta = ELECTRONS.get("mixing_beta", 0.4)

        smearing = t.cast(
            SmearingConfigurationSettingsModel,
            self.get_model("smearing"),
        )
        with smearing.hold_trait_notifications():
            if "degauss" in SYSTEM:
                smearing.degauss = SYSTEM["degauss"]
            if "smearing" in SYSTEM:
                smearing.type = SYSTEM["smearing"]

        num_kinds = len(self.input_structure.kinds) if self.has_structure else 1

        pseudos = t.cast(
            PseudosConfigurationSettingsModel,
            self.get_model("pseudos"),
        )
        with pseudos.hold_trait_notifications():
            if pseudo_family_string := state.get("pseudo_family"):
                pseudo_family = PseudoFamily.from_string(pseudo_family_string)
                library = f"{pseudo_family.library} {pseudo_family.accuracy}"
                if relativistic := pseudo_family.relativistic:
                    library += f" ({relativistic})"
                pseudos.functional = pseudo_family.functional
                pseudos.library = library
                pseudos.family = pseudo_family_string
            else:
                try:
                    pp_uuid = next(iter(PW["pseudos"].values()))
                    pseudo_info = get_pseudo_info(pp_uuid)
                    pseudos.functional = pseudo_info["functional"]
                except Exception:
                    pseudos.functional = None
                pseudos.library = None
                pseudos.family = None
                pseudos.show_upload_warning = True

            pseudos.functionals = [pseudos.functional] * num_kinds

            if pseudos_dictionary := PW.get("pseudos"):
                pseudos.dictionary = pseudos_dictionary
                pseudos.ecutwfc = SYSTEM.get("ecutwfc", 0.0)
                pseudos.ecutrho = SYSTEM.get("ecutrho", 0.0)

        magnetization = t.cast(
            MagnetizationConfigurationSettingsModel,
            self.get_model("magnetization"),
        )
        with magnetization.hold_trait_notifications():
            if "tot_magnetization" in SYSTEM:
                magnetization.type = "tot_magnetization"
                magnetization.total = SYSTEM["tot_magnetization"]

            if magnetic_moments := state.get("initial_magnetic_moments"):
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

        hubbard = t.cast(
            HubbardConfigurationSettingsModel,
            self.get_model("hubbard"),
        )
        with hubbard.hold_trait_notifications():
            if hubbard_parameters := state.get("hubbard_parameters"):
                hubbard.is_active = True
                hubbard.parameters = hubbard_parameters["hubbard_u"]
                starting_ns_eigenvalue = SYSTEM.get("starting_ns_eigenvalue")
                if starting_ns_eigenvalue is not None:
                    hubbard.set_active_eigenvalues(starting_ns_eigenvalue)

    def _link_model(self, model: PanelModel):
        tl.link(
            (self, "confirmed"),
            (model, "confirmed"),
        )
        super()._link_model(model)
