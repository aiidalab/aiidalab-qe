import ipywidgets as ipw

from aiidalab_qe.common.panel import ResultsPanel

from .model import WorkChainSummaryResultsModel


class WorkChainSummaryResultsPanel(ResultsPanel[WorkChainSummaryResultsModel]):
    title = "Workflow Summary"
    identifier = "summary"

    def render(self):
        if self.rendered:
            return
        report_html = self._model.generate_report_html()
        self.children = [ipw.HTML(report_html)]
        self.rendered = True
