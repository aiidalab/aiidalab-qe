import typing as t

import traitlets as tl

from aiidalab_qe.app.configuration import ConfigurationStepModel
from aiidalab_qe.app.result import ResultsStepModel
from aiidalab_qe.app.structure import StructureStepModel
from aiidalab_qe.app.submission import SubmissionStepModel
from aiidalab_qe.common.mixins import HasModels
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.wizard import State, WizardStepModel


class WizardModel(Model, HasModels[WizardStepModel]):
    state = tl.Dict(None, allow_none=True)
    selected_index = tl.Int(None, allow_none=True)
    loading_process = tl.Bool(False)

    def load_from_state(self, state: dict):
        step_index = state.get("step", 0) - 1
        structure_state = state.get("structure_state")
        configuration_state = state.get("configuration_state")
        resources_state = state.get("resources_state")
        process_uuid = state.get("process_uuid")

        if step_index >= 0:
            self.loading_process = True
            structure_model = t.cast(
                StructureStepModel,
                self.get_model("structure"),
            )
            structure_model.set_model_state(structure_state)

            if step_index >= 1:
                structure_model.confirm()
                configuration_model = t.cast(
                    ConfigurationStepModel,
                    self.get_model("configure"),
                )
                configuration_model.set_model_state(configuration_state)

                if step_index >= 2:
                    configuration_model.confirm()
                    submission_model = t.cast(
                        SubmissionStepModel,
                        self.get_model("submit"),
                    )
                    submission_model.set_model_state(resources_state)

                    if step_index >= 3:
                        submission_model.process_uuid = process_uuid
                        submission_model.confirm()
                        structure_model.lock()
                        configuration_model.lock()
                        submission_model.lock()

            self.loading_process = False

            self.selected_index = step_index

    def update_configuration_model(self):
        structure_model = t.cast(
            StructureStepModel,
            self.get_model("structure"),
        )
        configuration_model = t.cast(
            ConfigurationStepModel,
            self.get_model("configure"),
        )
        if structure_model.confirmed:
            configuration_model.await_properties()
            configuration_model.structure_uuid = structure_model.structure_uuid
        else:
            configuration_model.structure_uuid = None

    def update_submission_model(self):
        structure_model = t.cast(
            StructureStepModel,
            self.get_model("structure"),
        )
        configuration_model = t.cast(
            ConfigurationStepModel,
            self.get_model("configure"),
        )
        submission_model = t.cast(
            SubmissionStepModel,
            self.get_model("submit"),
        )
        if configuration_model.confirmed:
            submission_model.await_resources()
            submission_model.structure_uuid = structure_model.structure_uuid
            submission_model.input_parameters = configuration_model.get_model_state()
        else:
            submission_model.structure_uuid = None
            submission_model.input_parameters = {}

    def update_results_model(self):
        submission_model = t.cast(
            SubmissionStepModel,
            self.get_model("submit"),
        )
        results_model = t.cast(
            ResultsStepModel,
            self.get_model("results"),
        )
        tl.dlink(
            (submission_model, "process_uuid"),
            (results_model, "process_uuid"),
        )

    def lock_app(self):
        for identifier in ("structure", "configure", "submit"):
            model = self.get_model(identifier)
            model.lock()

    def auto_advance(self):
        if self.selected_index is None:
            return

        model_list = list(self._models.values())
        index = t.cast(int, self.selected_index)

        if (
            (selected_step := model_list[index]).auto_advance
            and not (index + 1 == len(model_list))
            and selected_step.state is State.SUCCESS
        ):
            self.selected_index += 1
