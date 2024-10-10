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
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)
from aiidalab_widgets_base import WizardAppWidgetStep

from .code import CodeModel, PluginCodes
from .model import SubmissionModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    def __init__(self, model: SubmissionModel, qe_auto_setup=True, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workflow submission panel")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_process_change,
            "process",
        )
        self._model.observe(
            self._on_input_parameters_change,
            "input_parameters",
        )
        self._model.observe(
            self._on_submission_blockers_change,
            [
                "internal_submission_blockers",
                "external_submission_blockers",
            ],
        )

        # # TODO for testing only - remove in PR
        # self._model.observe(
        #     lambda change: print(change["new"]),
        #     "input_parameters",
        # )

        self.qe_auto_setup = qe_auto_setup

        self._ALERT_MESSAGE = """
            <div class="alert alert-{alert_class} alert-dismissible">
                <a
                    href="#"
                    class="close"
                    data-dismiss="alert"
                    aria-label="close"
                >
                    &times;
                </a>
                <span
                    class="closebtn"
                    onclick="this.parentElement.style.display='none';"
                >
                    &times;
                </span>
                <strong>{message}</strong>
            </div>
        """

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.code_widgets_container = ipw.VBox()

        self.code_widgets: dict[str, QEAppComputationalResourcesWidget] = {}

        plugin_codes: PluginCodes = get_entry_items("aiidalab_qe.properties", "code")
        plugin_codes.update(
            {
                "dft": {
                    "pw": CodeModel(
                        description="pw.x:",
                        default_calc_job_plugin="quantumespresso.pw",
                        setup_widget_class=PwCodeResourceSetupWidget,
                    ),
                },
            }
        )
        for identifier, codes in plugin_codes.items():
            for name, code in codes.items():
                self._model.add_code(identifier, name, code)
                code.observe(
                    self._on_code_activation_change,
                    "is_active",
                )
                code.observe(
                    self._on_code_selection_change,
                    "selected",
                )

        plugin_codes["dft"]["pw"].activate()

        # set process label and description
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
        self.submit_button.on_click(self._on_submission)

        self.sssp_installation = PseudosInstallWidget()
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
        self.sssp_installation.observe(
            self._on_installation_change,
            ["busy", "installed"],
        )
        self.sssp_installation.observe(
            self._on_sssp_installed,
            "installed",
        )
        if self.qe_auto_setup:
            self.sssp_installation.refresh()

        self.qe_setup = QESetupWidget()
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
        self.qe_setup.observe(
            self._on_installation_change,
            ["busy", "installed"],
        )
        self.qe_setup.observe(
            self._on_qe_installed,
            "installed",
        )
        if self.qe_auto_setup:
            self.qe_setup.refresh()

        self.submission_blocker_messages = ipw.HTML()
        ipw.dlink(
            (self._model, "submission_blocker_messages"),
            (self.submission_blocker_messages, "value"),
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

    def set_submission_parameters(self, parameters):
        self._model.set_model_state(parameters)

    def set_process(self, process):
        self._model.process = process

    def reset(self):
        with self.hold_trait_notifications():
            self._model.reset()
            self._model.set_selected_codes()

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update()

    def _on_input_structure_change(self, _):
        self._update()

    def _on_process_change(self, _):
        with self.hold_trait_notifications():
            # TODO why here? Do we not populate traits earlier that would cover this?
            if self._model.process is not None:
                self._model.input_structure = self._model.process.inputs.structure
            self._update_state()

    def _on_input_parameters_change(self, _):
        self._update()

    def _on_submission_blockers_change(self, _):
        self._model.update_submission_blocker_message()
        self._update_state()

    def _on_installation_change(self, _):
        with self.hold_trait_notifications():
            self._model.update_submission_blockers()
            self._update_state()

    def _on_qe_installed(self, _):
        self._toggle_qe_installation_widget()

    def _on_sssp_installed(self, _):
        self._toggle_sssp_installation_widget()

    def _on_code_activation_change(self, change):
        self._toggle_code(change["owner"])

    def _on_code_selection_change(self, _):
        self._model.update_submission_blockers()

    def _on_submission(self, _):
        self._model.submit()
        self._update_state()

    def _toggle_sssp_installation_widget(self):
        sssp_installation_display = "none" if self._model.sssp_installed else "block"
        self.sssp_installation.layout.display = sssp_installation_display

    def _toggle_qe_installation_widget(self):
        qe_installation_display = "none" if self._model.qe_installed else "block"
        self.qe_setup.layout.display = qe_installation_display

    def _toggle_code(self, code: CodeModel):
        code_widget = code.get_setup_widget()
        code_widget.layout.display = "block" if code.is_active else "none"
        if not code.is_rendered:
            self._render_code_widget(code, code_widget)

    def _render_code_widget(self, code, code_widget):
        ipw.dlink(
            (code_widget.code_selection.code_select_dropdown, "options"),
            (code, "options"),
        )
        ipw.dlink(
            (code_widget, "value"),
            (code, "selected"),
            lambda value: value is not None,
        )
        code.observe(
            lambda change: setattr(code_widget, "parameters", change["new"]),
            "parameters",
        )
        self.code_widgets_container.children += (code_widget,)
        self._model.code_widgets[code.name] = code_widget
        self._model.set_selected_codes()
        code.is_rendered = True

    def _update(self):
        with self.hold_trait_notifications():
            self._model.update_active_codes()
            self._model.update_process_label()
            self._model.update_submission_blockers()

    def _update_state(self, _=None):
        if self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
        elif self.previous_step_state is not self.State.SUCCESS:
            self.state = self.State.INIT
        elif self._model.process is not None:
            self.state = self.State.SUCCESS
        else:
            if self._model.is_blocked:
                self.state = self.State.READY
            else:
                self._model.internal_submission_blockers = []
                self.state = self.state.CONFIGURED
