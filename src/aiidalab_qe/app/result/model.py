from __future__ import annotations

import contextlib
import time
import typing as t

import traitlets as tl

from aiida import orm
from aiida.engine import ProcessState
from aiida.engine.daemon.client import get_daemon_client
from aiida.engine.processes import control
from aiidalab_qe.common.mixins import HasProcess
from aiidalab_qe.common.process import STATE_ICONS
from aiidalab_qe.common.wizard import DependentWizardStepModel, State

from .utils import HasProcessModels, ResultsSubModel


class ResultsStepModel(
    DependentWizardStepModel,
    HasProcessModels[ResultsSubModel],
    HasProcess,
):
    identifier = "results"

    process_info = tl.Unicode("")
    process_remote_folder_is_clean = tl.Bool(False)
    kill_pending = tl.Bool(False)

    STATUS_TEMPLATE = "<h4>Workflow status: {}</h4"

    _dependencies = [
        "process_uuid",
    ]

    @property
    def is_active(self):
        return self.state is State.ACTIVE

    @property
    def is_failed(self):
        return self.state is State.FAIL

    def update_daemon_status(self):
        try:
            self.daemon_is_running = get_daemon_client().is_daemon_running
        except Exception:
            self.daemon_is_running = False
        self.daemon_status_known = True

    def kill_process(self):
        from aiidalab_qe.app.result.components.status.paused import PausedProcessesModel

        if not self.has_process:
            return

        paused_model = t.cast(PausedProcessesModel, self.get_model("status.paused"))
        if not self.daemon_status_known or not self.daemon_is_running:
            paused_model.error_message = (
                "The AiiDA daemon status is not available."
                if not self.daemon_status_known
                else "The AiiDA daemon is not running."
            )
            return

        if paused_model.paused_processes:
            paused_model.play_all()
            if paused_model.error_message:
                return
            self.kill_pending = True
            self._kill_deadline = time.monotonic() + 5.0
            return

        self._kill_root_workflow(paused_model)

    def process_pending_kill(self):
        from aiidalab_qe.app.result.components.status.paused import PausedProcessesModel

        if not self.kill_pending:
            return

        paused_model = t.cast(PausedProcessesModel, self.get_model("status.paused"))
        if not self.daemon_status_known or not self.daemon_is_running:
            paused_model.error_message = "The AiiDA daemon is not running."
            self.kill_pending = False
            return
        if paused_model.paused_processes:
            if time.monotonic() < self._kill_deadline:
                return
            paused_model.error_message = (
                "Could not resume all paused processes before killing the workflow."
            )
            self.kill_pending = False
            return

        self._kill_root_workflow(paused_model)

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
        self.daemon_is_running = False
        self.daemon_status_known = False
        self.kill_pending = False

    def _update(self, specific=""):
        self.update_daemon_status()
        self._update_process_remote_folder_state()

    def _update_process_remote_folder_state(self):
        if not (self.has_process and self.process.called_descendants):
            return
        cleaned = []
        for called_descendant in self.process.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    cleaned.append(called_descendant.outputs.remote_folder.is_empty)
        self.process_remote_folder_is_clean = all(cleaned)

    def _kill_root_workflow(self, paused_model):
        try:
            control.kill_processes([self.process])
        except Exception as exception:
            paused_model.error_message = str(exception)
            self.kill_pending = False
            return

        paused_model.reset()
        self.kill_pending = False

    def _get_process_status(self, state: str):
        return f"{state.capitalize()} {STATE_ICONS[state]}"
