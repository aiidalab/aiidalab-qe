from __future__ import annotations

from copy import deepcopy

import traitlets as tl
from aiida_pseudo.common.units import U

from aiida import orm
from aiida.common import exceptions
from aiida.plugins import GroupFactory
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.common.panel import PanelModel
from aiidalab_qe.setup.pseudos import PSEUDODOJO_VERSION, SSSP_VERSION, PseudoFamily

from .utils import get_pseudo_family_by_label, get_upf_dict

SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class PseudosConfigurationSettingsModel(
    PanelModel,
    HasInputStructure,
):
    title = "Pseudopotentials"
    identifier = "pseudos"

    dependencies = [
        "structure_uuid",
        "protocol",
        "spin_orbit",
    ]

    protocol = tl.Unicode()
    spin_orbit = tl.Unicode()

    # We allow `None` for the uuid in case the family does not contain a pseudo for a
    # given element, in which case, we block the app and notify the user.
    dictionary = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
        value_trait=tl.Unicode(allow_none=True),  # pseudopotential node uuid
        default_value={},
    )
    functional = tl.Unicode(allow_none=True)
    functional_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "PBE",
            "PBEsol",
        ],
    )
    functionals = tl.List(trait=tl.Unicode(allow_none=True))
    library = tl.Unicode(allow_none=True)
    library_options = tl.List(
        trait=tl.Unicode(),
        default_value=[
            "SSSP efficiency",
            "SSSP precision",
            "PseudoDojo standard (SR)",
            "PseudoDojo stringent (SR)",
            "PseudoDojo standard (FR)",
            "PseudoDojo stringent (FR)",
        ],
    )
    family = tl.Unicode(allow_none=True)
    family_header = tl.Unicode(allow_none=True)
    cutoffs = tl.List(
        trait=tl.List(tl.Float()),  # [[ecutwfc values], [ecutrho values]]
        default_value=[[0.0], [0.0]],
    )
    ecutwfc = tl.Float()
    ecutrho = tl.Float()
    status_message = tl.Unicode("", allow_none=True)
    show_upload_warning = tl.Bool(False)

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
            SSSP (v1.3) provides standard solid-state pseudopotentials, while PseudoDojo
            (v0.4) provides both scalar (SR) and full (FR) relativistic types.
            <br>
            If you are unsure, select 'SSSP efficiency', which for most calculations
            will produce sufficiently accurate results at comparatively small
            computational costs.
            <br>
            For higher accuracy, select 'SSSP precision' or 'PseudoDojo stringent'
            (computationally more expensive).
        </div>
    """

    family_help_message = tl.Unicode(PSEUDO_HELP_WO_SOC)

    include = True  # build-in panel

    def update(self, specific=""):
        if not self.has_structure:
            return
        family = self.family
        self.update_family_parameters()
        # When the app starts, the family is not yet set. `update_family_parameters`
        # will set the family, which will set the functionals and dictionary.
        # However, when the structure is changed, the family may already be set to
        # the default, in which case, the functionals and dictionary will not be
        # updated. Therefore, we need to force the update.
        if specific == "structure" and self.family == family:
            self.update_dictionary()
            self.update_functionals()

    def update_family_parameters(self):
        if self.locked:
            return

        if self.spin_orbit == "soc":
            if self.protocol in ["fast", "balanced"]:
                pseudo_family_string = (
                    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/FR/standard/upf"
                )
            else:
                pseudo_family_string = (
                    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/FR/stringent/upf"
                )
        else:
            protocol_inputs = PwBaseWorkChain.get_protocol_inputs(self.protocol)
            pseudo_family_string = protocol_inputs["pseudo_family"]

        pseudo_family = PseudoFamily.from_string(pseudo_family_string)
        library = f"{pseudo_family.library} {pseudo_family.accuracy}"
        if relativistic := pseudo_family.relativistic:
            library += f" ({relativistic})"

        self._defaults["functional"] = pseudo_family.functional
        self._defaults["library"] = library

        with self.hold_trait_notifications():
            self.functional = self._defaults["functional"]
            self.library = self._defaults["library"]

    def update_functionals(self):
        if self.locked or not (self.functional and self.family):
            return
        self.functionals = (
            [self.functional for _ in self.input_structure.kinds]
            if self.has_structure
            else []
        )

    def update_family(self):
        if self.locked or not (self.library and self.functional):
            return

        parts = self.library.split()
        library, accuracy = parts[:2]
        if len(parts) == 3:
            relativistic = parts[2].strip("()")
        functional = self.functional
        # XXX (jusong.yu): a validator is needed to check the family string is
        # consistent with the list of pseudo families defined in the setup_pseudos.py
        if library == "PseudoDojo":
            pseudo_family_string = f"PseudoDojo/{PSEUDODOJO_VERSION}/{functional}/{relativistic}/{accuracy}/upf"
        elif library == "SSSP":
            pseudo_family_string = f"SSSP/{SSSP_VERSION}/{functional}/{accuracy}"
        else:
            raise ValueError(
                f"Unknown pseudo family parameters: {library} | {accuracy}"
            )

        self._defaults["family"] = pseudo_family_string
        self.family = self._defaults["family"]

    def update_family_header(self):
        if not self.library:
            self.family_header = "<h4>Pseudopotential family</h4>"
            return

        library, accuracy = self.library.split()[:2]
        if library == "SSSP":
            pseudo_family_link = (
                f"https://www.materialscloud.org/discover/sssp/table/{accuracy}"
            )
        else:
            pseudo_family_link = "http://www.pseudo-dojo.org/"

        self.family_header = f"""
            <h4>
                <a href="{pseudo_family_link}" target="_blank">
                    Pseudopotential family
                </a>
            </h4>
        """

    def update_dictionary(self):
        if self.locked or not self.family:
            return

        if not self.has_structure:
            self._defaults |= {
                "dictionary": {},
                "cutoffs": [[0.0], [0.0]],
            }
            self.dictionary = self._get_default_dictionary()
            return

        pseudos = {}
        self.status_message = ""

        try:
            pseudo_family = get_pseudo_family_by_label(self.family)
        except Exception as err:
            raise ValueError(
                f"Failed to fetch pseudo family using the '{self.family}' string"
            ) from err

        pseudos = {}
        for kind in self.input_structure.kinds:
            # If the kind is not in the family, we set it to None.
            # This will block the app and notify the user of the missing pseudo.
            try:
                pseudo = pseudo_family.get_pseudo(kind.symbol)
            except ValueError:
                pseudo = None
            pseudos[kind.name] = pseudo

        self._defaults["dictionary"] = {
            kind: pseudo.uuid if pseudo else None for kind, pseudo in pseudos.items()
        }
        # Some pseudos may exist in more than one family (e.g., efficiency, precision).
        # This means that a change in the dictionary may not trigger an event due to
        # identical UUIDs. Therefore, we always reset the dictionary before updating it.
        # Note that dictionary observers bail if the dictionary is empty to avoid
        # unnecessary updates.
        self.dictionary = {}
        self.dictionary = self._get_default_dictionary()

    def update_cutoffs(self):
        if self.locked or not self.dictionary:
            return

        kinds = self.input_structure.kinds if self.has_structure else []
        self.status_message = ""

        if self.family:
            try:
                pseudo_family = get_pseudo_family_by_label(self.family)
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

        ecutwfc_list = []
        ecutrho_list = []
        for kind in kinds:
            if self.family:
                cutoff = cutoff_dict.get(kind.symbol, {})
                cutoff = {
                    key: U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                    for key, v in cutoff.items()
                }
                ecutwfc_list.append(cutoff.get("cutoff_wfc", 0.0))
                ecutrho_list.append(cutoff.get("cutoff_rho", 0.0))
            elif pp_uuid := self.dictionary.get(kind.name, None):
                try:
                    pp_node = orm.load_node(pp_uuid)
                    upf_dict = get_upf_dict(pp_node)
                except exceptions.NotExistent as err:
                    raise ValueError(
                        f"Pseudo potential with UUID {pp_uuid} does not exist"
                    ) from err
                except Exception as err:
                    raise ValueError(
                        f"Failed to read UPF data from pseudo potential with UUID {pp_uuid}"
                    ) from err

                try:
                    ecutwfc = float(int(upf_dict["header"]["wfc_cutoff"]))
                    ecutwfc_list.append(ecutwfc)
                except Exception:
                    ecutwfc_list.append(0.0)

                try:
                    ecutrho = float(int(upf_dict["header"]["rho_cutoff"]))
                    ecutrho_list.append(ecutrho)
                except KeyError:
                    ecutrho_list.append(0.0)

        self._defaults["cutoffs"] = [ecutwfc_list or [0.0], ecutrho_list or [0.0]]
        self.cutoffs = self._get_default_cutoffs()

    def update_library_options(self):
        if self.locked or not self.has_structure:
            return

        relativistic_options = [
            "PseudoDojo standard (FR)",
            "PseudoDojo stringent (FR)",
        ]
        if self.spin_orbit == "soc":
            library_options = relativistic_options
            self.family_help_message = self.PSEUDO_HELP_SOC
        else:
            library_options = [
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard (SR)",
                "PseudoDojo stringent (SR)",
                *relativistic_options,
            ]
            self.family_help_message = self.PSEUDO_HELP_WO_SOC

        self._defaults["library_options"] = library_options
        self.library_options = self._defaults["library_options"]

        self.update_family_parameters()

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
                    "PseudoDojo standard (SR)",
                    "PseudoDojo stringent (SR)",
                    "PseudoDojo standard (FR)",
                    "PseudoDojo stringent (FR)",
                ],
            )
        return super()._get_default(trait)

    def _get_default_dictionary(self):
        return deepcopy(self._defaults["dictionary"])

    def _get_default_cutoffs(self):
        return deepcopy(self._defaults["cutoffs"])

    def _check_blockers(self):
        if not (self.dictionary and self.functionals):
            return

        pseudos = []
        for kind_name, uuid in self.dictionary.items():
            kind = self.input_structure.get_kind(kind_name)
            try:
                assert uuid is not None
                pseudo = orm.load_node(uuid)
                pseudos.append(pseudo)
            except AssertionError:
                yield f"The selected pseudopotential family does not contain a pseudopotential for {kind.symbol}. Consider changing the family or uploading a custom pseudopotential."
                return
            except exceptions.NotExistent:
                yield f"Pseudopotential with UUID {uuid} does not exist for {kind.symbol}."
                return

        functional_set = set(self.functionals)
        if len(functional_set) != 1:
            yield "All pseudopotentials must have the same exchange-correlation (XC) functional."
        elif self.functional and self.functional not in functional_set:
            yield "Selected exchange-correlation (XC) functional is not consistent with the pseudopotentials."

        relativistic_set = {pp.base.extras.get("relativistic", None) for pp in pseudos}
        if self.spin_orbit == "soc":
            if relativistic_set != {"full"}:
                yield "For spin-orbit coupling (SOC) calculations, all pseudopotentials must be fully relativistic."
