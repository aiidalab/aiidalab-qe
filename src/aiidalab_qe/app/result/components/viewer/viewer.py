from __future__ import annotations

import asyncio
import html
from collections.abc import Callable
from importlib.metadata import distributions

import ipywidgets as ipw
from packaging.utils import canonicalize_name

from aiidalab_qe.app.result.components import ResultsComponent
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ResultsPanel
from aiidalab_qe.plugins.registry import (
    get_plugin_registry,
    get_plugin_statuses,
)
from aiidalab_qe.plugins.utils import get_entry_items

from .model import WorkflowResultsViewerModel
from .structure import StructureResultsModel, StructureResultsPanel


class WorkflowResultsViewer(ResultsComponent[WorkflowResultsViewerModel]):
    def __init__(
        self,
        model: WorkflowResultsViewerModel,
        check_plugin_compatibility: bool = True,
        on_plugin_status_checked: Callable[[set[str]], None] | None = None,
        **kwargs,
    ):
        self.check_plugin_compatibility = check_plugin_compatibility
        self.on_plugin_status_checked = on_plugin_status_checked
        self.plugin_status_checked = False
        # NOTE: here we want to add the structure and plugin models to the viewer
        # model BEFORE we define the observation of the process uuid. This ensures
        # that when the process changes, its reflected in the sub-models prior to
        # the logic of the process change event handler.
        # TODO avoid exceptions! Ensure sub-model synchronization in general!
        self.panels: dict[str, ResultsPanel] = {}
        self.incompatible_plugin_results: list[str] = []
        self.incompatible_plugin_ids: set[str] = set()
        self._add_structure_panel(model)
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
        self.plugin_compatibility_warning = ipw.HTML(
            value=(
                self._get_plugin_compatibility_warning()
                if self.plugin_status_checked
                else "Loading/validating plugin results..."
            )
        )
        self.children = [
            InAppGuide(identifier="results-panel"),
            self.plugin_compatibility_warning,
            self.tabs,
        ]
        self.ensure_plugin_status_check()

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

    def ensure_plugin_status_check(self) -> None:
        """Start runtime plugin checks once the viewer is on screen."""
        if self.plugin_status_checked or getattr(self, "_plugin_status_task", None):
            return
        if not self.check_plugin_compatibility:
            self.plugin_status_checked = True
            self.plugin_compatibility_warning.value = ""
            self._notify_plugin_status_checked()
            return
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            asyncio.run(self._fetch_plugin_results_async())
        else:
            self._plugin_status_task = loop.create_task(
                self._fetch_plugin_results_async()
            )

    async def _fetch_plugin_results_async(self) -> None:
        registry = get_plugin_registry()
        statuses = await get_plugin_statuses(registry.plugins)
        self._fetch_plugin_results(self._model, registry, statuses)
        self.plugin_status_checked = True
        if self.rendered:
            self.plugin_compatibility_warning.value = (
                self._get_plugin_compatibility_warning()
            )
            self._set_tabs()
        self._notify_plugin_status_checked()

    def _notify_plugin_status_checked(self) -> None:
        if self.on_plugin_status_checked is not None:
            self.on_plugin_status_checked(self.incompatible_plugin_ids)

    def _fetch_plugin_results(self, viewer_model, registry, statuses):
        incompatible_packages = {}
        for plugin_name, plugin_data in registry.plugins.items():
            status = statuses[plugin_name]
            if status.is_incompatible:
                incompatible_packages[status.distribution_name] = plugin_data.title

        for distribution in distributions():
            plugin_name = distribution.metadata.get("Name")
            if not plugin_name or canonicalize_name(plugin_name) not in (
                incompatible_packages
            ):
                continue

            self.incompatible_plugin_ids.update(
                entry_point.name
                for entry_point in distribution.entry_points
                if entry_point.group == "aiidalab_qe.properties"
            )

        def is_entry_compatible(entry_point):
            distribution = getattr(entry_point, "dist", None)
            metadata = getattr(distribution, "metadata", {})
            package_name = metadata.get("Name") if metadata else None
            title = (
                incompatible_packages.get(canonicalize_name(package_name))
                if package_name
                else None
            )
            if title:
                self.incompatible_plugin_results.append(title)
                return False
            return True

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
            panel = entry["panel"]
            model = entry["model"]()
            viewer_model.add_model(identifier, model)
            self.panels[identifier] = panel(model=model)

    def _get_plugin_compatibility_warning(self) -> str:
        plugins = sorted(set(self.incompatible_plugin_results))
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
