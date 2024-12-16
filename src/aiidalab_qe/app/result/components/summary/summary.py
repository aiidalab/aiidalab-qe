import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponent

from .model import WorkChainSummaryModel
from .outputs import WorkChainOutputs


class WorkChainSummary(ResultsComponent[WorkChainSummaryModel]):
    def __init__(self, model: WorkChainSummaryModel, **kwargs):
        super().__init__(model=model, **kwargs)
        self.has_report = False
        self.has_output = False

    def _render(self):
        self._render_summary()

    def _on_process_change(self, _):
        if not self.has_report:
            self._render_summary()

    def _on_monitor_counter_change(self, _):
        if not self.has_output:
            self._render_output()

    def _render_summary(self):
        if not self._model.has_process:
            return
        report = self._model.generate_report_html()
        self.children = [ipw.HTML(report)]
        self.has_report = True

    def _render_output(self):
        process_node = self._model.fetch_process_node()
        if process_node and process_node.is_terminated:
            self.children += (WorkChainOutputs(node=process_node),)
            self.has_output = True
