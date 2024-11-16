import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel

from .model import WorkChainSummaryModel


class WorkChainSummary(ResultsPanel[WorkChainSummaryModel]):
    title = "Workflow Summary"
    identifier = "summary"

    def render(self):
        if self.rendered:
            return
        report_html = self._model.generate_report_html()
        self.children = [ipw.HTML(report_html)]
        self.rendered = True
