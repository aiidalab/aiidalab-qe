from __future__ import annotations

from copy import deepcopy

import traitlets as tl
from aiida_pseudo.common.units import U

from aiida.common import exceptions
from aiida.plugins import GroupFactory
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION, PseudoFamily
from aiidalab_qe.utils import fetch_pseudo_family_by_label

from ..subsettings import AdvancedCalculationSubSettingsModel

SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class PseudosConfigurationSettingsModel(
    AdvancedCalculationSubSettingsModel,
    HasInputStructure,
):
    identifier = "pseudos"

    dependencies = [
        "input_structure",
        "protocol",
        "spin_orbit",
    ]

    protocol = tl.Unicode()
    spin_orbit = tl.Unicode()

    dictionary = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
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
    functional = tl.Unicode(DEFAULT["advanced"]["pseudo_family"]["functional"])
    functional_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "PBE",
            "PBEsol",
        ],
    )
    library = tl.Unicode(
        " ".join(
            [
                DEFAULT["advanced"]["pseudo_family"]["library"],
                DEFAULT["advanced"]["pseudo_family"]["accuracy"],
            ]
        )
    )
    library_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "SSSP efficiency",
            "SSSP precision",
            "PseudoDojo standard",
            "PseudoDojo stringent",
        ],
    )
    cutoffs = tl.List(
        trait=tl.List(tl.Float()),  # [[ecutwfc values], [ecutrho values]]
        default_value=[[0.0], [0.0]],
    )
    ecutwfc = tl.Float()
    ecutrho = tl.Float()
    status_message = tl.Unicode("")

    PSEUDO_HELP_SOC = """
        <div class="pseudo-text">
            Spin-orbit coupling (SOC) calculations are supported exclusively with
            PseudoDojo pseudopotentials. PseudoDojo offers these pseudopotentials
            in two versions: standard and stringent. Here, we utilize the FR
            (fully relativistic) type from PseudoDojo. Please ensure you choose
            appropriate cutoff values for your calculations.
        </div>
    """

    PSEUDO_HELP_WO_SOC = """
        <div class="pseudo-text">
            If you are unsure, select 'SSSP efficiency', which for most calculations
            will produce sufficiently accurate results at comparatively small
            computational costs.
            <br>
            If your calculations require a higher accuracy, select 'SSSP accuracy' or
            'PseudoDojo stringent', which will be computationally more expensive.
            <br>
            SSSP is the standard solid-state pseudopotentials.
            The PseudoDojo version used here is the SR relativistic type.
        </div>
    """

    family_help_message = tl.Unicode(PSEUDO_HELP_WO_SOC)

    pseudo_filename_reset_trigger = tl.Int(0)

    def update(self, specific=""):  # noqa: ARG002
        with self.hold_trait_notifications():
            if not self.has_structure:
                self._defaults |= {
                    "dictionary": {},
                    "cutoffs": [[0.0], [0.0]],
                }
            else:
                self.update_default_pseudos()
                self.update_default_cutoffs()
            self.update_family_parameters()
            self.update_family()

    def update_default_pseudos(self):
        if self.loaded_from_process:
            return

        pseudos = {}
        self.status_message = ""

        try:
            pseudo_family = fetch_pseudo_family_by_label(self.family)
            pseudos = pseudo_family.get_pseudos(structure=self.input_structure)
        except ValueError as exception:
            self.status_message = f"""
                <div class='alert alert-danger'>
                    ERROR: {exception!s}
                </div>
            """
            return

        self._defaults["dictionary"] = {
            kind: pseudo.uuid for kind, pseudo in pseudos.items()
        }
        self.dictionary = self._get_default_dictionary()

    def update_default_cutoffs(self):
        """Update wavefunction and density cutoffs from pseudo family."""
        if self.loaded_from_process:
            return

        kinds = []
        self.status_message = ""

        try:
            pseudo_family = fetch_pseudo_family_by_label(self.family)
            current_unit = pseudo_family.get_cutoffs_unit()
            cutoff_dict = pseudo_family.get_cutoffs()
        except exceptions.NotExistent:
            self.status_message = f"""
                <div class='alert alert-danger'>
                    ERROR: required pseudo family `{self.family}` is
                    not installed. Please use `aiida-pseudo install` to install
                    it."
                </div>
            """
        except ValueError as exception:
            self.status_message = f"""
                <div class='alert alert-danger'>
                    ERROR: failed to obtain recommended cutoffs for pseudos
                    `{pseudo_family}`: {exception}
                </div>
            """
        else:
            kinds = self.input_structure.kinds if self.input_structure else []

        ecutwfc_list = []
        ecutrho_list = []
        for kind in kinds:
            cutoff = cutoff_dict.get(kind.symbol, {})
            cutoff = {
                key: U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                for key, v in cutoff.items()
            }
            ecutwfc_list.append(cutoff["cutoff_wfc"])
            ecutrho_list.append(cutoff["cutoff_rho"])

        self._defaults["cutoffs"] = [ecutwfc_list or [0.0], ecutrho_list or [0.0]]
        self.cutoffs = self._get_default_cutoffs()

    def update_library_options(self):
        if self.loaded_from_process:
            return
        if self.spin_orbit == "soc":
            library_options = [
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.family_help_message = self.PSEUDO_HELP_SOC
        else:
            library_options = [
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.family_help_message = self.PSEUDO_HELP_WO_SOC
        self._defaults["library_options"] = library_options
        self.library_options = self._defaults["library_options"]

        self.update_family_parameters()

    def update_family_parameters(self):
        if self.loaded_from_process:
            return
        if self.spin_orbit == "soc":
            if self.protocol in ["fast", "moderate"]:
                pseudo_family_string = "PseudoDojo/0.4/PBEsol/FR/standard/upf"
            else:
                pseudo_family_string = "PseudoDojo/0.4/PBEsol/FR/stringent/upf"
        else:
            pseudo_family_string = PwBaseWorkChain.get_protocol_inputs(self.protocol)[
                "pseudo_family"
            ]

        pseudo_family = PseudoFamily.from_string(pseudo_family_string)

        self._defaults["library"] = f"{pseudo_family.library} {pseudo_family.accuracy}"
        self._defaults["functional"] = pseudo_family.functional

        with self.hold_trait_notifications():
            self.library = self._defaults["library"]
            self.functional = self._defaults["functional"]

    def update_family(self):
        if self.loaded_from_process:
            return
        library, accuracy = self.library.split()
        functional = self.functional
        # XXX (jusong.yu): a validator is needed to check the family string is
        # consistent with the list of pseudo families defined in the setup_pseudos.py
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

        self._defaults["family"] = pseudo_family_string
        self.family = self._defaults["family"]

    def reset(self):
        with self.hold_trait_notifications():
            self.dictionary = self._get_default("dictionary")
            self.cutoffs = self._get_default("cutoffs")
            self.ecutwfc = max(self.cutoffs[0])
            self.ecutrho = max(self.cutoffs[1])
            self.library_options = self._get_default("library_options")
            self.library = self._get_default("library")
            self.functional = self._get_default("functional")
            self.functional_options = self._get_default("functional_options")
            self.family = self._get_default("family")
            self.family_help_message = self._get_default("family_help_message")
            self.status_message = self._get_default("status_message")
            self.pseudo_filename_reset_trigger += 1

    def _get_default(self, trait):
        if trait == "dictionary":
            return deepcopy(self._defaults.get(trait, {}))
        if trait == "cutoffs":
            return deepcopy(self._defaults.get(trait, [[0.0], [0.0]]))
        if trait == "functional_options":
            return self._defaults.get(
                trait,
                [
                    "PBE",
                    "PBEsol",
                ],
            )
        if trait == "library_options":
            return self._defaults.get(
                trait,
                [
                    "SSSP efficiency",
                    "SSSP precision",
                    "PseudoDojo standard",
                    "PseudoDojo stringent",
                ],
            )
        return self._defaults.get(trait, self.traits()[trait].default_value)

    def _get_default_dictionary(self):
        return deepcopy(self._defaults["dictionary"])

    def _get_default_cutoffs(self):
        return deepcopy(self._defaults["cutoffs"])
