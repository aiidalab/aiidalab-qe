from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiida_quantumespresso.common.types import RelaxType
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import (
    HasInputStructure,
    HasModels,
)
from aiidalab_qe.common.panel import ConfigurationSettingsModel
from aiidalab_qe.common.wizard import QeConfirmableWizardStepModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore

NO_RELAXATION_OPTION = ("Structure as is", "none")


class ConfigurationStepModel(
    QeConfirmableWizardStepModel,
    HasModels[ConfigurationSettingsModel],
    HasInputStructure,
):
    identifier = "configuration"

    relax_type_help = tl.Unicode()
    relax_type_options = tl.List([NO_RELAXATION_OPTION])
    relax_type = tl.Unicode(NO_RELAXATION_OPTION[-1], allow_none=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._default_models = {
            "workchain",
            "advanced",
        }

        self.relax_type_help_template = """
            <div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
                You have {option_count} options:
                <br>
                (1) Structure as is: perform a self consistent calculation using
                the structure provided as input.
                <br>
                (2) Atomic positions: perform a full relaxation of the internal
                atomic coordinates.
                {full_relaxation_option}
            </div>
        """

    def update(self):
        if self.has_pbc:
            relax_type_help = self.relax_type_help_template.format(
                option_count="three",
                full_relaxation_option=(
                    """
                    <br>
                    (3) Full geometry: perform a full relaxation of the internal atomic
                    coordinates and the cell parameters.
                    """
                ),
            )
            relax_type_options = [
                NO_RELAXATION_OPTION,
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ]
        else:
            relax_type_help = self.relax_type_help_template.format(
                option_count="two",
                full_relaxation_option="",
            )
            relax_type_options = [
                NO_RELAXATION_OPTION,
                ("Atomic positions", "positions"),
            ]

        default = DEFAULT["workchain"]["relax_type"]
        default_available = [default in [option[1] for option in relax_type_options]]
        relax_type = default if default_available else relax_type_options[-1][-1]

        self._defaults = {
            "relax_type_help": relax_type_help,
            "relax_type_options": relax_type_options,
            "relax_type": relax_type,
        }

        self.relax_type_help = self._get_default_relax_type_help()
        self.relax_type_options = self._get_default_relax_type_options()
        self.relax_type = self._get_default_relax_type()

    def get_model_state(self):
        parameters = {
            identifier: model.get_model_state()
            for identifier, model in self.get_models()
            if model.include
        }
        parameters["workchain"] |= {
            "relax_type": self.relax_type,
            "properties": self._get_properties(),
        }
        return parameters

    def set_model_state(self, parameters):
        with self.hold_trait_notifications():
            workchain_parameters: dict = parameters.get("workchain", {})
            self.relax_type = workchain_parameters.get("relax_type")
            properties = set(workchain_parameters.get("properties", []))
            for identifier, model in self.get_models():
                model.include = identifier in self._default_models | properties
                if parameters.get(identifier):
                    model.set_model_state(parameters[identifier])
                    model.loaded_from_process = True

    def reset(self):
        self.confirmed = False
        self.relax_type_help = self._get_default_relax_type_help()
        self.relax_type_options = self._get_default_relax_type_options()
        self.relax_type = self._get_default_relax_type()
        for identifier, model in self.get_models():
            if identifier not in self._default_models:
                model.include = False

    def _link_model(self, model: ConfigurationSettingsModel):
        ipw.link(
            (self, "confirmed"),
            (model, "confirmed"),
        )
        super()._link_model(model)

    def _get_properties(self):
        properties = []
        for identifier, model in self.get_models():
            if identifier in self._default_models:
                continue
            if model.include:
                properties.append(identifier)
        if RelaxType(self.relax_type) is not RelaxType.NONE or not properties:
            properties.append("relax")
        return properties

    def _get_default_relax_type_help(self):
        return self._defaults.get("relax_type_help", "")

    def _get_default_relax_type_options(self):
        return self._defaults.get("relax_type_options", [NO_RELAXATION_OPTION])

    def _get_default_relax_type(self):
        options = self._get_default_relax_type_options()
        relax_type = self._defaults.get("relax_type", NO_RELAXATION_OPTION[-1])
        return (
            relax_type
            if relax_type in [option[1] for option in options]
            else options[-1][-1]
        )
