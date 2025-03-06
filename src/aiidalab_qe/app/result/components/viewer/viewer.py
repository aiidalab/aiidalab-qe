from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ResultsPanel

from .model import WorkChainResultsViewerModel
from .structure import StructureResultsModel, StructureResultsPanel


class WorkChainResultsViewer(ResultsComponent[WorkChainResultsViewerModel]):
    def __init__(self, model: WorkChainResultsViewerModel, **kwargs):
        # NOTE: here we want to add the structure and plugin models to the viewer
        # model BEFORE we define the observation of the process uuid. This ensures
        # that when the process changes, its reflected in the sub-models prior to
        # the logic of the process change event handler.
        # TODO avoid exceptions! Ensure sub-model synchronization in general!
        self.panels: dict[str, ResultsPanel] = {}
        self._add_structure_panel(model)
        self._fetch_plugin_results(model)
        super().__init__(model=model, **kwargs)

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
        self.tabs = ipw.Tab(selected_index=None)
        self.tabs.observe(
            self._on_tab_change,
            "selected_index",
        )
        self.children = [
            InAppGuide(identifier="results-panel"),
            self.tabs,
        ]

    def _post_render(self):
        self._set_tabs()

    def _update_panels(self):
        self.panels = {
            identifier: self.panels[identifier]
            for identifier, model in self._model.get_models()
            if model.include
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

    def _add_structure_panel(self, viewer_model: WorkChainResultsViewerModel):
        structure_model = StructureResultsModel()
        structure_model.process_uuid = viewer_model.process_uuid
        self.structure_results = StructureResultsPanel(model=structure_model)
        identifier = structure_model.identifier
        viewer_model.add_model(identifier, structure_model)
        self.panels = {
            identifier: self.structure_results,
            **self.panels,
        }

    def _fetch_plugin_results(self, viewer_model: WorkChainResultsViewerModel):
        entries = get_entry_items("aiidalab_qe.properties", "result")
        for identifier, entry in entries.items():
            for key in ("panel", "model"):
                if key not in entry:
                    raise ValueError(
                        f"Entry {identifier} is missing the results '{key}' key"
                    )
            panel = entry["panel"]
            model = entry["model"]()
            viewer_model.add_model(identifier, model)
            self.panels[identifier] = panel(
                identifier=identifier,
                model=model,
            )
