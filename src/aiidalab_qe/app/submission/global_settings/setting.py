"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw

from aiidalab_qe.common.code import CodeModel, PluginCodes
from aiidalab_qe.common.panel import ResourceSettingsPanel
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget

from .model import GlobalResourceSettingsModel


class GlobalResourceSettingsPanel(ResourceSettingsPanel[GlobalResourceSettingsModel]):
    def __init__(self, model: GlobalResourceSettingsModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_input_parameters_change,
            "input_parameters",
        )
        self._model.observe(
            self._on_plugin_overrides_change,
            "plugin_overrides",
        )

    def render(self):
        if self.rendered:
            return

        self.code_widgets_container = ipw.VBox()

        self.plugin_overrides_notification = ipw.HTML()
        ipw.dlink(
            (self._model, "plugin_overrides_notification"),
            (self.plugin_overrides_notification, "value"),
        )

        self.children = [
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
                    Use the other tabs if you want to override the resource settings per plugin.
                </div>
            """),
            self.code_widgets_container,
            self.plugin_overrides_notification,
        ]

        self.rendered = True

        # Render any active codes
        for _, code_model in self._model.get_models():
            if code_model.is_active:
                self._toggle_code(code_model)

    def set_up_codes(self, codes: PluginCodes):
        for identifier, code_models in codes.items():
            for _, code_model in code_models.items():
                base_code_model = self._model.add_global_model(identifier, code_model)
                if base_code_model is not None:
                    base_code_model.observe(
                        self._on_code_activation_change,
                        "is_active",
                    )
                    base_code_model.observe(
                        self._on_code_selection_change,
                        "selected",
                    )
        self._model.update_global_codes()

    def reset(self):
        self._model.set_selected_codes()

    def _on_input_parameters_change(self, _):
        self._model.update_active_codes()

    def _on_input_structure_change(self, _):
        self._model.check_resources()

    def _on_plugin_overrides_change(self, _):
        self._model.update_plugin_overrides_notification()

    def _on_code_activation_change(self, change):
        self._toggle_code(change["owner"])

    def _on_code_selection_change(self, _):
        self._model.update_blockers()

    def _on_pw_code_resource_change(self, _):
        self._model.check_resources()

    def _on_code_resource_change(self, _):
        self._model.update_global_codes()

    def _render_code_widget(
        self,
        code_model: CodeModel,
        code_widget: QEAppComputationalResourcesWidget,
    ):
        super()._render_code_widget(code_model, code_widget)
        if code_model.default_calc_job_plugin == "quantumespresso.pw":
            code_model.observe(
                self._on_pw_code_resource_change,
                [
                    "num_cpus",
                    "num_nodes",
                    "ntasks_per_node",
                    "cpus_per_task",
                    "max_wallclock_seconds",
                ],
            )

        def toggle_widget(_=None, model=code_model, widget=code_widget):
            widget = self.code_widgets[model.name]
            widget.layout.display = "block" if model.is_active else "none"

        code_model.observe(
            toggle_widget,
            "is_active",
        )

        toggle_widget()
