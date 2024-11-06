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
from aiidalab_widgets_base import WizardAppWidgetStep

from .code import CodeModel, PluginCodes, PwCodeModel
from .model import SubmissionModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    def __init__(self, model: SubmissionModel, qe_auto_setup=True, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading workflow submission panel")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_submission,
            "confirmed",
        )
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_input_parameters_change,
            "input_parameters",
        )
        self._model.observe(
            self._on_process_change,
            "process",
        )
        self._model.observe(
            self._on_submission_blockers_change,
            [
                "internal_submission_blockers",
                "external_submission_blockers",
            ],
        )
        self._model.observe(
            self._on_installation_change,
            ["installing_sssp", "sssp_installed"],
        )
        self._model.observe(
            self._on_sssp_installed,
            "sssp_installed",
        )
        self._model.observe(
            self._on_installation_change,
            ["installing_qe", "qe_installed"],
        )
        self._model.observe(
            self._on_qe_installed,
            "qe_installed",
        )

        self.code_widgets: dict[str, QEAppComputationalResourcesWidget] = {}

        self._install_sssp(qe_auto_setup)
        self._set_up_qe(qe_auto_setup)
        self._set_up_codes()

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.code_widgets_container = ipw.VBox()

        self.process_label = ipw.Text(
            description="Label:",
            layout=ipw.Layout(width="auto", indent="0px"),
        )
        ipw.link(
            (self._model, "process_label"),
            (self.process_label, "value"),
        )
        self.process_description = ipw.Textarea(
            description="Description",
            layout=ipw.Layout(width="auto", indent="0px"),
        )
        ipw.link(
            (self._model, "process_description"),
            (self.process_description, "value"),
        )

        self.submit_button = ipw.Button(
            description="Submit",
            tooltip="Submit the calculation with the selected parameters.",
            icon="play",
            button_style="success",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )
        ipw.dlink(
            (self, "state"),
            (self.submit_button, "disabled"),
            lambda state: state != self.State.CONFIGURED,
        )
        self.submit_button.on_click(self.submit)

        self.submission_blocker_messages = ipw.HTML()
        ipw.dlink(
            (self._model, "submission_blocker_messages"),
            (self.submission_blocker_messages, "value"),
        )

        self.submission_warning_messages = ipw.HTML()
        ipw.dlink(
            (self._model, "submission_warning_messages"),
            (self.submission_warning_messages, "value"),
        )

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Codes</h4>
                </div>
            """),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
                    Select the code to use for running the calculations. The codes on
                    the local machine (localhost) are installed by default, but you can
                    configure new ones on potentially more powerful machines by clicking
                    on "Setup new code".
                </div>
            """),
            self.code_widgets_container,
            self.sssp_installation,
            self.qe_setup,
            self.submission_blocker_messages,
            self.submission_warning_messages,
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Labeling Your Job</h4>
                    <p style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
                        Label your job and provide a brief description. These details
                        help identify the job later and make the search process easier.
                        While optional, adding a description is recommended for better
                        clarity.
                    </p>
                </div>
            """),
            self.process_label,
            self.process_description,
            self.submit_button,
        ]

        self.rendered = True

        # Render any active codes
        self._model.get_code("dft", "pw").activate()
        for _, code in self._model.get_code_models(flat=True):
            if code.is_active:
                self._toggle_code(code)

    def submit(self, _=None):
        self._model.confirm()

    def reset(self):
        with self.hold_trait_notifications():
            self._model.reset()
            self._model.set_selected_codes()

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update_state()

    def _on_input_structure_change(self, _):
        self._model.check_resources()

    def _on_input_parameters_change(self, _):
        self._model.update_active_codes()
        self._model.update_process_label()
        self._model.update_submission_blockers()

    def _on_process_change(self, _):
        with self.hold_trait_notifications():
            # TODO why here? Do we not populate traits earlier that would cover this?
            if self._model.process_node is not None:
                self._model.input_structure = self._model.process_node.inputs.structure
            self._update_state()

    def _on_submission_blockers_change(self, _):
        self._model.update_submission_blocker_message()
        self._update_state()

    def _on_installation_change(self, _):
        self._model.update_submission_blockers()

    def _on_qe_installed(self, _):
        self._toggle_qe_installation_widget()
        if self._model.qe_installed:
            self._model.refresh_codes()

    def _on_sssp_installed(self, _):
        self._toggle_sssp_installation_widget()

    def _on_code_activation_change(self, change):
        self._toggle_code(change["owner"])

    def _on_code_selection_change(self, _):
        self._model.update_submission_blockers()

    def _on_pw_code_resource_change(self, _):
        self._model.check_resources()

    def _on_submission(self, _):
        self._update_state()

    def _install_sssp(self, qe_auto_setup):
        self.sssp_installation = PseudosInstallWidget(auto_start=False)
        ipw.dlink(
            (self.sssp_installation, "busy"),
            (self._model, "installing_sssp"),
        )
        ipw.dlink(
            (self.sssp_installation, "installed"),
            (self._model, "installing_sssp"),
            lambda installed: not installed,
        )
        ipw.dlink(
            (self.sssp_installation, "installed"),
            (self._model, "sssp_installed"),
        )
        if qe_auto_setup:
            self.sssp_installation.refresh()

    def _set_up_qe(self, qe_auto_setup):
        self.qe_setup = QESetupWidget(auto_start=False)
        ipw.dlink(
            (self.qe_setup, "busy"),
            (self._model, "installing_qe"),
        )
        ipw.dlink(
            (self.qe_setup, "installed"),
            (self._model, "installing_qe"),
            lambda installed: not installed,
        )
        ipw.dlink(
            (self.qe_setup, "installed"),
            (self._model, "qe_installed"),
        )
        if qe_auto_setup:
            self.qe_setup.refresh()

    def _set_up_codes(self):
        codes: PluginCodes = {
            "dft": {
                "pw": PwCodeModel(
                    description="pw.x",
                    default_calc_job_plugin="quantumespresso.pw",
                    code_widget_class=PwCodeResourceSetupWidget,
                ),
            },
            **get_entry_items("aiidalab_qe.properties", "code"),
        }
        for identifier, code_models in codes.items():
            for name, code_model in code_models.items():
                self._model.add_code(identifier, name, code_model)
                code_model.observe(
                    self._on_code_activation_change,
                    "is_active",
                )
                code_model.observe(
                    self._on_code_selection_change,
                    "selected",
                )

    def _toggle_sssp_installation_widget(self):
        sssp_installation_display = "none" if self._model.sssp_installed else "block"
        self.sssp_installation.layout.display = sssp_installation_display

    def _toggle_qe_installation_widget(self):
        qe_installation_display = "none" if self._model.qe_installed else "block"
        self.qe_setup.layout.display = qe_installation_display

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
        ipw.dlink(
            (code_model, "options"),
            (code_widget.code_selection.code_select_dropdown, "options"),
        )
        ipw.link(
            (code_model, "selected"),
            (code_widget.code_selection.code_select_dropdown, "value"),
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
        if isinstance(code_widget, PwCodeResourceSetupWidget):
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
                ["num_cpus", "num_nodes"],
            )
        code_widgets = self.code_widgets_container.children[:-1]  # type: ignore
        self.code_widgets_container.children = [*code_widgets, code_widget]
        code_model.is_rendered = True

    def _update_state(self, _=None):
        if self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
        elif self.previous_step_state is not self.State.SUCCESS:
            self.state = self.State.INIT
        elif self._model.process_node is not None:
            self.state = self.State.SUCCESS
        elif self._model.is_blocked:
            self.state = self.State.READY
        else:
            self.state = self.state.CONFIGURED
