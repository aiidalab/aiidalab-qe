import typing as t

import ipywidgets as ipw

from aiidalab_qe.common.mixins import HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import LoadingWidget


class ResultsComponentModel(Model, HasProcess):
    identifier = "results"


RCM = t.TypeVar("RCM", bound=ResultsComponentModel)


class ResultsComponent(ipw.VBox, t.Generic[RCM]):
    def __init__(self, model: RCM, **kwargs):
        self.loading_message = LoadingWidget(f"Loading {model.identifier}")

        super().__init__(
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_process_change,
            "process_uuid",
        )
        self._model.observe(
            self._on_monitor_counter_change,
            "monitor_counter",
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return
        self._render()
        self.rendered = True
        self._post_render()

    def _on_process_change(self, _):
        pass

    def _on_monitor_counter_change(self, _):
        pass

    def _render(self):
        raise NotImplementedError

    def _post_render(self):
        pass
