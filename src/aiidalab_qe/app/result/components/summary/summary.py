import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponent

from .model import WorkChainSummaryModel
from .outputs import WorkChainOutputs


class WorkChainSummary(ResultsComponent[WorkChainSummaryModel]):
    def __init__(self, model: WorkChainSummaryModel, **kwargs):
        super().__init__(model=model, **kwargs)
        self.has_report = False
        self.has_output = False

    def _on_process_change(self, _):
        if not self.has_report:
            self._render_summary()

    def _on_monitor_counter_change(self, _):
        if not self.has_output:
            self._render_output()

    def _render(self):
        self._render_summary()

    def _render_summary(self):
        if not self._model.has_process:
            return

        settings_summary = ipw.HTML(
            value=self._model.generate_report_html(),
        )
        settings_summary.add_class("summary-panel")

        self.output_download_container = ipw.VBox(
            children=[
                ipw.HTML("""
                    <div style="line-height: 140%; margin: 0; margin-bottom: 10px;">
                        <h2>Download the data</h2>
                        Once the workflow is finished, you can download raw data
                        (i.e. input and output files) and/or the AiiDA archive
                        (ready to be shared or imported into another AiiDA profile).
                    </div>
                """),
                ipw.HTML("Download buttons will appear here when available."),
            ],
        )
        self.output_download_container.add_class("summary-panel")

        self.container = ipw.HBox(
            children=[
                settings_summary,
                self.output_download_container,
            ],
        )
        self.container.add_class("workflow-summary-container")
        self.children = [self.container]
        self.has_report = True

    def _render_output(self):
        process_node = self._model.fetch_process_node()
        if process_node and process_node.is_terminated:
            output_download_widget = WorkChainOutputs(node=process_node)
            output_download_widget.layout.width = "100%"
            self.output_download_container.children = [
                self.output_download_container.children[0],  # type: ignore
                output_download_widget,
            ]
            self.has_output = True
