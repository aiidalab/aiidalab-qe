import typing as t

import ipywidgets as ipw

from aiidalab_qe.app.result.utils import ResultsSubModel
from aiidalab_widgets_base import LoadingWidget

RSM = t.TypeVar("RSM", bound=ResultsSubModel)


class ResultsComponent(ipw.VBox, t.Generic[RSM]):
    def __init__(self, model: RSM, **kwargs):
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
