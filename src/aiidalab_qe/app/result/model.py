from __future__ import annotations

import contextlib

import traitlets as tl

from aiida import orm
from aiida.engine import ProcessState
from aiida.engine.processes import control
from aiidalab_qe.common.mixins import HasModels, HasProcess
from aiidalab_qe.common.process import STATE_ICONS
from aiidalab_qe.common.wizard import DependentWizardStepModel, State

from .components import ResultsComponentModel


class ResultsStepModel(
    DependentWizardStepModel,
    HasModels[ResultsComponentModel],
    HasProcess,
):
    identifier = "results"

    process_info = tl.Unicode("")
    process_remote_folder_is_clean = tl.Bool(False)

    STATUS_TEMPLATE = "<h4>Workflow status: {}</h4"

    _dependencies = [
        "process_uuid",
    ]

    @property
    def is_loading(self):
        return super().is_loading and self.state is not State.ACTIVE

    def update(self):
        self._update_process_remote_folder_state()

    def kill_process(self):
        if self.has_process:
            control.kill_processes([self.process])

    def clean_remote_data(self):
        if not self.has_process:
            return
        for called_descendant in self.process.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    called_descendant.outputs.remote_folder._clean()
        self.process_remote_folder_is_clean = True

    def update_state(self):
        if not self.has_process:
            self.state = State.INIT
            return

        if process_state := self.process.process_state:
            status = self._get_process_status(process_state.value)
        else:
            status = "Unknown"

        if process_state is ProcessState.CREATED:
            self.state = State.ACTIVE
        elif process_state in (
            ProcessState.RUNNING,
            ProcessState.WAITING,
        ):
            self.state = State.ACTIVE
            status = self._get_process_status("running")  # overwrite status
        elif process_state in (
            ProcessState.EXCEPTED,
            ProcessState.KILLED,
        ):
            self.state = State.FAIL
        elif self.process.is_failed:
            self.state = State.FAIL
        elif self.process.is_finished_ok:
            self.state = State.SUCCESS
        else:
            self.state = State.CONFIGURED

        self.process_info = self.STATUS_TEMPLATE.format(status)

    def reset(self):
        self.process_uuid = None
        self.process_info = ""

    def _update_process_remote_folder_state(self):
        if not (self.has_process and self.process.called_descendants):
            return
        cleaned = []
        for called_descendant in self.process.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    cleaned.append(called_descendant.outputs.remote_folder.is_empty)
        self.process_remote_folder_is_clean = all(cleaned)

    def _link_model(self, model: ResultsComponentModel):
        tl.dlink(
            (self, "process_uuid"),
            (model, "process_uuid"),
        )
        tl.dlink(
            (self, "monitor_counter"),
            (model, "monitor_counter"),
        )

    def _get_process_status(self, state: str):
        return f"{state.capitalize()} {STATE_ICONS[state]}"
