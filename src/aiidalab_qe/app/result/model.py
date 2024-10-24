from __future__ import annotations

import contextlib

import traitlets as tl

from aiida import orm
from aiida.engine.processes import control


class ResultsModel(tl.HasTraits):
    process = tl.Unicode(allow_none=True)

    process_info = tl.Unicode("")
    process_remote_folder_is_clean = tl.Bool(False)

    _process_node: orm.ProcessNode | None = None

    @property
    def process_node(self):
        if self._process_node is None:
            self._process_node = self._get_process_node()
        return self._process_node

    def update(self):
        self._update_process_remote_folder_state()

    def kill_process(self):
        if process := self._get_process_node():
            control.kill_processes(process)

    def clean_remote_data(self):
        if self.process_node is None:
            return
        for called_descendant in self.process_node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    called_descendant.outputs.remote_folder._clean()
        self.process_remote_folder_is_clean = True

    def reset(self):
        self.process = None
        self.process_info = ""

    def _get_process_node(self):
        return orm.load_node(self.process) if self.process else None

    def _update_process_remote_folder_state(self):
        if self.process_node is None:
            return
        cleaned = []
        for called_descendant in self.process_node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                with contextlib.suppress(Exception):
                    cleaned.append(called_descendant.outputs.remote_folder.is_empty)
        self.process_remote_folder_is_clean = all(cleaned)
