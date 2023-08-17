import ipywidgets as ipw

from .report import generate_report_html


class SummaryView(ipw.VBox):
    def __init__(self, wc_node, **kwargs):
        report_html = generate_report_html(wc_node)

        self.summary_view = ipw.HTML(report_html)
        super().__init__(
            children=[self.summary_view],
            **kwargs,
        )
