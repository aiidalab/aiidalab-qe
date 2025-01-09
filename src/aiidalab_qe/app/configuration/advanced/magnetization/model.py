from copy import deepcopy

import traitlets as tl
from aiida_pseudo.groups.family import PseudoPotentialFamily

from aiida.common.exceptions import NotExistent
from aiida_quantumespresso.workflows.protocols.utils import get_starting_magnetization
from aiidalab_qe.common.mixins import HasInputStructure
from aiidalab_qe.utils import fetch_pseudo_family_by_label

from ..subsettings import AdvancedCalculationSubSettingsModel


class MagnetizationConfigurationSettingsModel(
    AdvancedCalculationSubSettingsModel,
    HasInputStructure,
):
    identifier = "magnetization"

    dependencies = [
        "input_structure",
        "electronic_type",
        "spin_type",
        "pseudos.family",
    ]

    electronic_type = tl.Unicode()
    spin_type = tl.Unicode()
    family = tl.Unicode()

    type_options = tl.List(
        trait=tl.List(tl.Unicode()),
        default_value=[
            ["Initial magnetic moments", "starting_magnetization"],
            ["Total magnetization", "tot_magnetization"],
        ],
    )
    type = tl.Unicode("starting_magnetization")
    type_help = tl.Unicode("")
    total = tl.Float(0.0)
    moments = tl.Dict(
        key_trait=tl.Unicode(),  # kind name
        value_trait=tl.Float(),  # magnetic moment
        default_value={},
    )

    _TYPE_HELP_TEMPLATE = """
        <div style="line-height: 1.4; margin-bottom: 5px">
            {content}
        </div>
    """

    def update(self, specific=""):  # noqa: ARG002
        if self.spin_type == "none" or not self.has_structure:
            self._defaults["moments"] = {}
        else:
            self.update_type_help()
            self._update_default_moments()
        with self.hold_trait_notifications():
            self.moments = self._get_default_moments()

    def update_type_help(self):
        if self.electronic_type == "insulator" or self.type == "tot_magnetization":
            self.type_help = self._TYPE_HELP_TEMPLATE.format(
                content="""
                    Constrain the desired total electronic magnetization (difference
                    between majority and minority spin charge).
                """
            )
        else:
            self.type_help = self._TYPE_HELP_TEMPLATE.format(
                content="""
                    If a nonzero ground-state magnetization is expected, you
                    <strong>must</strong> assign a nonzero value to at least one atomic
                    type (note that the app might already provide tentative initial
                    values to chemical elements that typically display magnetism).
                    To simulate an antiferromagnetic state, first, if you have not
                    done so already, please use the atom tag editor <b>(Select
                    structure -> Edit structure -> Edit atom tags)</b> to mark atoms of
                    the species of interest as distinct by assigning each a different
                    integer tag. Once tagged, assign each an initial magnetic moments
                    of opposite sign.
                """
            )

    def reset(self):
        with self.hold_trait_notifications():
            self.type = self.traits()["type"].default_value
            self.total = self.traits()["total"].default_value
            self.moments = self._get_default_moments()

    def _update_default_moments(self):
        try:
            family = fetch_pseudo_family_by_label(self.family)
            guess = get_starting_magnetization(self.input_structure, family)
            self._defaults["moments"] = {
                kind.name: self._to_moment(kind.symbol, guess[kind.name], family)
                for kind in self.input_structure.kinds
            }
        except NotExistent:
            self._defaults["moments"] = {
                kind_name: 0.0 for kind_name in self.input_structure.get_kind_names()
            }

    @staticmethod
    def _to_moment(
        symbol: str,
        magnetization: float,
        family: PseudoPotentialFamily,
    ):
        return round(magnetization * family.get_pseudo(symbol).z_valence, 3)

    def _get_default_moments(self):
        return deepcopy(self._defaults.get("moments", {}))
