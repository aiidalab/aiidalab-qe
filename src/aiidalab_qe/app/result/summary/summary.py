import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel

from .model import WorkChainSummaryModel


class WorkChainSummary(ResultsPanel[WorkChainSummaryModel]):
    title = "Workflow Summary"
    identifier = "summary"

    def render(self):
        if self.rendered:
            return
        self._update_view()
        self.rendered = True

    def _update_view(self):
        self.report_html = self._model.generate_report_html()
        self.children = [
            ipw.HTML(self.report_html),
        ]

    def _on_monitor_counter_change(self, _):
        pass

    def _update_process_state_notification(self):
        pass
