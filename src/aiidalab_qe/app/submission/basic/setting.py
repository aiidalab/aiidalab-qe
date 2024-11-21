"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.common.setup_pseudos import PseudosInstallWidget
from aiidalab_qe.common.widgets import (
    LoadingWidget,
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)
from aiidalab_qe.common.code import CodeModel, PluginCodes, PwCodeModel
from aiidalab_qe.common.panel import SettingsPanel
from .model import BasicCodeModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class BasicCodeSettings(SettingsPanel[BasicCodeModel]):
    title = "Basic"
    identifier = "basic"


    def __init__(self, model: BasicCodeModel, **kwargs):
        super().__init__(model, **kwargs)
        self._set_up_codes()
        self._model.observe(
            self._on_input_parameters_change,
            "input_parameters",
        )
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

    def render(self):
        if self.rendered:
            return

        self.code_widgets_container = ipw.VBox()
        self.code_widgets = {}

        self.children = [
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
                    Use the other tabs if you want to override the resource settings per plugin.
                </div>
            """),
            self.code_widgets_container,
        ]

        self.rendered = True
        # Render any active codes
        self._model.get_code("quantumespresso.pw").activate()
        for code_model in self._model.codes.values():
            if code_model.is_active:
                self._toggle_code(code_model)

    def reset(self):
        with self.hold_trait_notifications():
            self._model.reset()
            self._model.set_selected_codes()

    def _on_input_parameters_change(self, _):
        self._model.update_active_codes()
    
    def _on_input_structure_change(self, _):
        self._model.check_resources()

    def _on_code_activation_change(self, change):
        self._toggle_code(change["owner"])

    def _on_code_selection_change(self, _):
        """"""
        #TODO: update the selected code in the input parameters
        # self._model.update_submission_blockers()

    def _on_pw_code_resource_change(self, _):
        self._model.check_resources()

    def _on_code_resource_change(self, _):
        """Update the plugin resources."""
        # trigger the update of basic codes
        self._model.basic_codes = {}
        self._model.basic_codes = self._model.get_model_state()["codes"]

    def _set_up_codes(self):
        
        codes: PluginCodes = {
            "dft": {
                "pw": PwCodeModel(
                    description="pw.x",
                    default_calc_job_plugin="quantumespresso.pw",
                    code_widget_class=PwCodeResourceSetupWidget,
                ),
            },
        }
        # Load codes from plugins
        eps = get_entry_items("aiidalab_qe.properties", "code")
        for identifier, data in eps.items():
            codes[identifier] = data["model"].codes
        for identifier, code_models in codes.items():
            for _, code_model in code_models.items():
                # use the new code model created using the basic code model
                code_model = self._model.add_code(identifier, code_model)
                if code_model is not None:
                    code_model.observe(
                        self._on_code_activation_change,
                        "is_active",
                    )
                    code_model.observe(
                        self._on_code_selection_change,
                        "selected",
                    )

    def _toggle_code(self, code_model: CodeModel):
        if not self.rendered:
            return
        if not code_model.is_rendered:
            loading_message = LoadingWidget(f"Loading {code_model.name} code")
            self.code_widgets_container.children += (loading_message,)
        if code_model.name not in self.code_widgets:
            code_widget = code_model.code_widget_class(
                description=code_model.description,
                default_calc_job_plugin=code_model.default_calc_job_plugin,
            )
            self.code_widgets[code_model.name] = code_widget
        else:
            code_widget = self.code_widgets[code_model.name]
        code_widget.layout.display = "block" if code_model.is_active else "none"
        if not code_model.is_rendered:
            self._render_code_widget(code_model, code_widget)

    def _render_code_widget(
        self,
        code_model: CodeModel,
        code_widget: QEAppComputationalResourcesWidget,
    ):
        code_model.update()
        ipw.dlink(
            (code_model, "options"),
            (code_widget.code_selection.code_select_dropdown, "options"),
        )
        ipw.link(
            (code_model, "selected"),
            (code_widget.code_selection.code_select_dropdown, "value"),
        )
        code_widget.code_selection.code_select_dropdown.observe(
            self._on_code_selection_change,
            "value",
        )
        ipw.dlink(
            (code_model, "selected"),
            (code_widget.code_selection.code_select_dropdown, "disabled"),
            lambda selected: not selected,
        )
        ipw.link(
            (code_model, "num_cpus"),
            (code_widget.num_cpus, "value"),
        )
        ipw.link(
            (code_model, "num_nodes"),
            (code_widget.num_nodes, "value"),
        )
        ipw.link(
            (code_model, "ntasks_per_node"),
            (code_widget.resource_detail.ntasks_per_node, "value"),
        )
        ipw.link(
            (code_model, "cpus_per_task"),
            (code_widget.resource_detail.cpus_per_task, "value"),
        )
        ipw.link(
            (code_model, "max_wallclock_seconds"),
            (code_widget.resource_detail.max_wallclock_seconds, "value"),
        )
        if code_model.default_calc_job_plugin == "quantumespresso.pw":
            ipw.link(
                (code_model, "override"),
                (code_widget.parallelization.override, "value"),
            )
            ipw.link(
                (code_model, "npool"),
                (code_widget.parallelization.npool, "value"),
            )
            code_model.observe(
                self._on_pw_code_resource_change,
                ["num_cpus", "num_nodes", "ntasks_per_node", "cpus_per_task", "max_wallclock_seconds"],
            )
        code_model.observe(
            self._on_code_resource_change,
            ["num_cpus", "num_nodes", "ntasks_per_node", "cpus_per_task", "max_wallclock_seconds"],
        )
        code_widgets = self.code_widgets_container.children[:-1]  # type: ignore
        self.code_widgets_container.children = [*code_widgets, code_widget]
        code_model.is_rendered = True
