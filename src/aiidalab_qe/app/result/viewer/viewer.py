from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.app.result.summary.model import WorkChainSummaryModel
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_widgets_base import register_viewer_widget

from ..structure import StructureResults, StructureResultsModel
from ..summary import WorkChainSummary
from .model import WorkChainViewerModel
from .outputs import WorkChainOutputs


@register_viewer_widget("process.workflow.workchain.WorkChainNode.")
class WorkChainViewer(ipw.VBox):
    _results_shown = tl.Set()

    def __init__(self, node, model: WorkChainViewerModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading result panels")],
            **kwargs,
        )

        self._model = model
        self._model.process_uuid = node.uuid

        self.rendered = False

        summary_model = WorkChainSummaryModel()
        summary_model.process_uuid = node.uuid
        self.summary = WorkChainSummary(model=summary_model)
        self._model.add_model("summary", summary_model)

        self.results: dict[str, ResultsPanel] = {
            "summary": self.summary,
        }

        # TODO consider refactoring structure relaxation panel as a plugin
        if "relax" in self._model.properties:
            self._add_structure_panel()

        self._fetch_plugin_results()

    def render(self):
        if self.rendered:
            return

        node = self._model.fetch_process_node()

        self.title = ipw.HTML()

        title = "<hr style='height:2px;background-color:#2097F3;'>"
        if node:
            formula = node.inputs.structure.get_formula()
            title += f"\n<h4>QE App Workflow (pk: {node.pk}) &mdash; {formula}</h4>"
        else:
            title += "\n<h4>QE App Workflow</h4>"

        self.title.value = title

        self.tabs = ipw.Tab(selected_index=None)

        self.children = [
            self.title,
            self.tabs,
        ]

        self.rendered = True

        self._update_tabs()

        if node and node.is_finished:
            self._add_workflow_output_widget()

    def _update_tabs(self):
        children = []
        titles = []
        for identifier, model in [*self._model.get_models()]:
            if model.include:
                results = self.results[identifier]
                titles.append(results.title)
                children.append(results)
        self.tabs.children = children
        for i, title in enumerate(titles):
            self.tabs.set_title(i, title)
        self.tabs.selected_index = 0
        self.summary.render()

    def _add_workflow_output_widget(self):
        process_node = self._model.fetch_process_node()
        self.summary.children += (WorkChainOutputs(node=process_node),)

    def _add_structure_panel(self):
        structure_model = StructureResultsModel()
        structure_model.process_uuid = self._model.process_uuid
        self.structure_results = StructureResults(model=structure_model)
        identifier = self.structure_results.identifier
        self._model.add_model(identifier, structure_model)
        self.results[identifier] = self.structure_results

    def _fetch_plugin_results(self):
        entries = get_entry_items("aiidalab_qe.properties", "result")
        for identifier, entry in entries.items():
            for key in ("panel", "model"):
                if key not in entry:
                    raise ValueError(
                        f"Entry {identifier} is missing the results '{key}' key"
                    )
            panel = entry["panel"]
            model = entry["model"]()
            model.process_uuid = self._model.process_uuid
            self.results[identifier] = panel(
                identifier=identifier,
                model=model,
            )
            self._model.add_model(identifier, model)
