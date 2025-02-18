from __future__ import annotations

import contextlib

import traitlets as tl

from aiida import orm
from aiida.engine.processes import control
from aiidalab_qe.common.mixins import HasModels, HasProcess
from aiidalab_qe.common.wizard import QeWizardStepModel

from .components import ResultsComponentModel


class ResultsStepModel(
    QeWizardStepModel,
    HasModels[ResultsComponentModel],
    HasProcess,
):
    identifier = "results"

    process_info = tl.Unicode("")
    process_remote_folder_is_clean = tl.Bool(False)

    def update(self):
        self._update_process_remote_folder_state()

    def kill_process(self):
        if process_node := self.fetch_process_node():
            control.kill_processes([process_node])

    def clean_remote_data(self):
        if not (process_node := self.fetch_process_node()):
            return
        for called_descendant in process_node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    called_descendant.outputs.remote_folder._clean()
        self.process_remote_folder_is_clean = True

    def reset(self):
        self.process_uuid = None
        self.process_info = ""

    def _update_process_remote_folder_state(self):
        process_node = self.fetch_process_node()
        if not (process_node and process_node.called_descendants):
            return
        cleaned = []
        for called_descendant in process_node.called_descendants:
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
