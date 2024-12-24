from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.panel import ResultsPanel

from .model import WorkChainResultsViewerModel
from .structure import StructureResults, StructureResultsModel


class WorkChainResultsViewer(ResultsComponent[WorkChainResultsViewerModel]):
    def __init__(self, model: WorkChainResultsViewerModel, **kwargs):
        super().__init__(model=model, **kwargs)
        self.panels: dict[str, ResultsPanel] = {}
        self._fetch_plugin_results()

    def _on_process_change(self, _):
        self._update_panels()
        if self.rendered:
            self._set_tabs()

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: ResultsPanel = self.tabs.children[tab_index]  # type: ignore
        tab.render()

    def _render(self):
        if node := self._model.fetch_process_node():
            formula = node.inputs.structure.get_formula()
            title = f"\n<h4>QE App Workflow (pk: {node.pk}) &mdash; {formula}</h4>"
        else:
            title = "\n<h4>QE App Workflow</h4>"

        self.title = ipw.HTML(title)

        self.tabs = ipw.Tab(selected_index=None)
        self.tabs.observe(
            self._on_tab_change,
            "selected_index",
        )

        # TODO consider refactoring structure relaxation panel as a plugin
        if "relax" in self._model.properties:
            self._add_structure_panel()

        self.children = [
            self.title,
            self.tabs,
        ]

    def _post_render(self):
        self._set_tabs()

    def _update_panels(self):
        properties = self._model.properties
        need_electronic_structure = "bands" in properties and "pdos" in properties
        self.panels = {
            identifier: panel
            for identifier, panel in self.panels.items()
            if identifier in properties
            or (identifier == "electronic_structure" and need_electronic_structure)
        }

    def _set_tabs(self):
        children = []
        titles = []
        for identifier, model in self._model.get_models():
            if identifier not in self.panels:
                continue
            results = self.panels[identifier]
            titles.append(model.title)
            children.append(results)
        self.tabs.children = children
        for i, title in enumerate(titles):
            self.tabs.set_title(i, title)
        if children:
            self.tabs.selected_index = 0

    def _add_structure_panel(self):
        structure_model = StructureResultsModel()
        structure_model.process_uuid = self._model.process_uuid
        self.structure_results = StructureResults(model=structure_model)
        identifier = structure_model.identifier
        self._model.add_model(identifier, structure_model)
        self.panels = {
            identifier: self.structure_results,
            **self.panels,
        }

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
            self._model.add_model(identifier, model)
            self.panels[identifier] = panel(
                identifier=identifier,
                model=model,
            )
