from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
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

    def __init__(self, node: orm.Node, model: WorkChainViewerModel, **kwargs):
        if node.process_label != "QeAppWorkChain":
            super().__init__()
            return

        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading result panels")],
            **kwargs,
        )

        self._model = model

        self.rendered = False

        summary_model = WorkChainSummaryModel()
        self.summary = WorkChainSummary(model=summary_model)
        self._model.add_model("summary", summary_model)

        self.results: dict[str, ResultsPanel] = {
            "summary": self.summary,
        }

        # TODO consider refactoring structure relaxation panel as a plugin
        if "relax" in self._model.properties:
            self._add_structure_panel()

        self._fetch_plugin_results()

        # HACK
        self.render()

    def render(self):
        if self.rendered:
            return

        node = self._model.process_node

        self.title = ipw.HTML(f"""
            <hr style="height:2px;background-color:#2097F3;">
            <h4>
                QE App Workflow (pk: {node.pk}) &mdash; {node.inputs.structure.get_formula()}
            </h4>
        """)

        self.tabs = ipw.Tab(selected_index=None)
        self.tabs.observe(
            self._on_tab_change,
            "selected_index",
        )

        self.children = [
            self.title,
            self.tabs,
        ]

        self.rendered = True

        self._update_tabs()

        if node.is_finished:
            self._add_workflow_output_widget()

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: ResultsPanel = self.tabs.children[tab_index]  # type: ignore
        tab.render()

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

    def _add_workflow_output_widget(self):
        self.summary.children += (WorkChainOutputs(self._model.process_node),)

    def _add_structure_panel(self):
        structure_model = StructureResultsModel()
        self.structure_results = StructureResults(model=structure_model)
        identifier = self.structure_results.identifier
        self._model.add_model(identifier, structure_model)
        self.results[identifier] = self.structure_results

    def _fetch_plugin_results(self):
        entries = get_entry_items("aiidalab_qe.properties", "result")
        needs_electronic_structure = all(
            identifier in self._model.properties for identifier in ("bands", "pdos")
        )
        for identifier, entry in entries.items():
            for key in ("panel", "model"):
                if key not in entry:
                    raise ValueError(
                        f"Entry {identifier} is missing the results '{key}' key"
                    )
            panel = entry["panel"]
            model = entry["model"]()
            self._model.add_model(identifier, model)
            if identifier == "electronic_structure" and needs_electronic_structure:
                model.include = True
            else:
                model.include = identifier in self._model.properties
            self.results[identifier] = panel(
                identifier=identifier,
                model=model,
            )
