import ipywidgets as ipw

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.result.components import ResultsComponent

from .model import WorkChainSummaryModel
from .outputs import WorkChainOutputs

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class WorkChainSummary(ResultsComponent[WorkChainSummaryModel]):
    has_settings_report = False
    has_download_widget = False

    def _on_process_change(self, _):
        if not self.has_settings_report:
            self._render_summary()

    def _on_monitor_counter_change(self, _):
        if not self.rendered:
            return
        if not self.has_download_widget:
            self._render_download_widget()
        if not self._model.has_failure_report:
            self._model.generate_failure_report()

    def _render(self):
        self._render_summary()

    def _render_summary(self):
        if not self._model.has_process:
            return

        settings_summary = ipw.HTML(
            value=self._model.generate_report_html(),
        )
        settings_summary.add_class("summary-panel")

        self.output_download_help = ipw.HTML("""
            <div style="line-height: 1.4; margin: 0; margin-bottom: 10px;">
                <h2>Download the data</h2>
                Once the workflow is finished, you can download raw data (i.e. input
                and output files) and/or the AiiDA archive (ready to be shared or
                imported into another AiiDA profile).
            </div>
        """)

        self.output_download_container = ipw.VBox(
            children=[
                self.output_download_help,
                ipw.HTML("""
                    <div style="line-height: 1.4">
                        Download buttons will appear here when available.
                    </div>
                """),
            ],
        )
        self.output_download_container.add_class("summary-panel")
        self._render_download_widget()

        container = ipw.HBox(
            children=[
                settings_summary,
                self.output_download_container,
            ],
        )
        container.add_class("workflow-summary-container")
        container.add_class(DEFAULT["summary_format"])

        self.failed_calculation_report = ipw.HTML()
        ipw.dlink(
            (self._model, "failed_calculation_report"),
            (self.failed_calculation_report, "value"),
        )

        self.children = [
            container,
            self.failed_calculation_report,
        ]
        self.has_settings_report = True

    def _render_download_widget(self):
        process_node = self._model.fetch_process_node()
        if process_node and process_node.is_terminated:
            output_download_widget = WorkChainOutputs(node=process_node)
            output_download_widget.layout.width = "100%"
            self.output_download_container.children = [
                self.output_download_help,
                output_download_widget,
            ]
            self.has_download_widget = True
