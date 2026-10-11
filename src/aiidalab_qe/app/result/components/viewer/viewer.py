from __future__ import annotations

import html

import ipywidgets as ipw

from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.plugins.utils import get_entry_items

from .model import WorkflowResultsViewerModel
from .structure import StructureResultsModel, StructureResultsPanel

_VALIDATION_MESSAGE = """
    <div class="loading" style="display: flex; align-items: center; font-size: unset;">
        Loading/validating plugins...
        <i class="fa fa-spinner fa-spin fa-1x fa-fw"></i>
    </div>
    """


class WorkflowResultsViewer(ResultsComponent[WorkflowResultsViewerModel]):
    def __init__(self, model: WorkflowResultsViewerModel, **kwargs):
        self.panels: dict[str, ResultsPanel] = {}
        self._result_plugins_loaded = False
        self._add_structure_panel(model)
        super().__init__(model=model, **kwargs)
        self._model.observe(
            self._on_plugin_status_checked,
            "plugin_status_checked",
        )

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
        self.plugin_compatibility_warning = ipw.HTML(
            value=(
                self._get_plugin_compatibility_warning()
                if self._model.plugin_status_checked
                else _VALIDATION_MESSAGE
            )
        )
        self.children = [
            InAppGuide(identifier="results-panel"),
            self.plugin_compatibility_warning,
            self.tabs,
        ]

    def _post_render(self):
        if self._model.plugin_status_checked:
            self._load_result_plugins()
        self._set_tabs()
        self.plugin_compatibility_warning.value = (
            self._get_plugin_compatibility_warning()
            if self._model.plugin_status_checked
            else _VALIDATION_MESSAGE
        )

    def _update_panels(self):
        self.panels = {
            identifier: self.panels[identifier]
            for identifier, model in self._model.get_models()
            if model.include and identifier in self.panels
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
        for index, title in enumerate(titles):
            self.tabs.set_title(index, title)
        if children:
            self.tabs.selected_index = 0

    def _add_structure_panel(self, viewer_model: WorkflowResultsViewerModel):
        structure_model = StructureResultsModel()
        structure_model.process_uuid = viewer_model.process_uuid
        self.structure_results = StructureResultsPanel(model=structure_model)
        identifier = structure_model.identifier
        viewer_model.add_model(identifier, structure_model)
        self.panels = {
            identifier: self.structure_results,
            **self.panels,
        }

    def _on_plugin_status_checked(self, change):
        if not change["new"]:
            return
        self._load_result_plugins()
        if self.rendered:
            self.plugin_compatibility_warning.value = (
                self._get_plugin_compatibility_warning()
            )
            self._set_tabs()

    def _load_result_plugins(self):
        if self._result_plugins_loaded:
            return

        incompatible_plugin_ids = self._model.incompatible_plugin_ids

        def is_entry_compatible(entry_point):
            return entry_point.name not in incompatible_plugin_ids

        entries = get_entry_items(
            "aiidalab_qe.properties",
            "result",
            entry_point_filter=is_entry_compatible,
        )
        for identifier, entry in entries.items():
            for key in ("panel", "model"):
                if key not in entry:
                    raise ValueError(
                        f"Entry {identifier} is missing the results '{key}' key"
                    )
            model = entry["model"]()
            if self._model.process_uuid:
                model.process_uuid = self._model.process_uuid
                model.update()
                if not model.include:
                    continue
            self._model.add_model(identifier, model)
            self.panels[identifier] = entry["panel"](model=model)
        self._result_plugins_loaded = True

    def _get_plugin_compatibility_warning(self) -> str:
        plugins = sorted(set(self._model.incompatible_plugin_titles))
        if not plugins:
            return ""
        items = "".join(f"<li>{html.escape(title)}</li>" for title in plugins)
        return (
            '<div class="alert alert-danger" role="alert">'
            "Result panels are disabled for the following plugins:"
            f"<ul>{items}</ul>"
            'Visit the <a href="./plugin_manager.ipynb" target="_blank">Plugin store</a> to '
            "update or remove these plugins, then reload the app.</div>"
        )
