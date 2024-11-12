from __future__ import annotations

import contextlib

import traitlets as tl

from aiida import orm
from aiida.engine.processes import control
from aiidalab_qe.common.mixins import HasProcess
from aiidalab_qe.common.mvc import Model


class ResultsStepModel(Model, HasProcess):
    process_info = tl.Unicode("")
    process_remote_folder_is_clean = tl.Bool(False)

    def update(self):
        self._update_process_remote_folder_state()

    def kill_process(self):
        if self.has_process:
            control.kill_processes([self.process_node])

    def clean_remote_data(self):
        if not self.has_process:
            return
        for called_descendant in self.process_node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    called_descendant.outputs.remote_folder._clean()
        self.process_remote_folder_is_clean = True

    def reset(self):
        self.process_uuid = None
        self.process_info = ""

    def _update_process_remote_folder_state(self):
        if not self.has_process:
            return
        cleaned = []
        for called_descendant in self.process_node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    cleaned.append(called_descendant.outputs.remote_folder.is_empty)
        self.process_remote_folder_is_clean = all(cleaned)
