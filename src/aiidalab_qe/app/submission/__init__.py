"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.engine import ProcessBuilderNamespace
from aiida.engine import submit as aiida_submit
from aiida.orm.utils.serialize import serialize
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.common.setup_pseudos import PseudosInstallWidget
from aiidalab_qe.common.widgets import (
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)
from aiidalab_qe.workflows import QeAppWorkChain
from aiidalab_widgets_base import WizardAppWidgetStep

from .code import CodeModel, PluginCodes
from .model import SubmissionModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    internal_submission_blockers = tl.List(tl.Unicode())
    external_submission_blockers = tl.List(tl.Unicode())

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

    @property
    def is_blocked(self):
        return any(
            [
                *self.internal_submission_blockers,
                *self.external_submission_blockers,
            ]
        )

    def render(self):
        if self.rendered:
            return

        self.message_area = ipw.Output()

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
            self.message_area,
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

    @tl.observe("internal_submission_blockers", "external_submission_blockers")
    def _on_submission_blockers_change(self, _):
        self._update_submission_blocker_message()
        self._update_state()

    def _on_input_structure_change(self, _):
        self._update()

    def _on_process_change(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            # TODO why here? Do we not populate traits earlier that would cover this?
            if process_node is not None:
                self._model.input_structure = process_node.inputs.structure
            self._update_state()

    def _on_input_parameters_change(self, _):
        self._update()

    def _on_submission(self, _):
        self._submit()
        self._update_state()

    def _on_installation_change(self, _):
        with self.hold_trait_notifications():
            self._update_submission_blockers()
            self._update_state()

    def _on_qe_installed(self, change):
        self._model.installing_qe = not change["new"]
        self.qe_setup.layout.display = "none" if change["new"] else "block"

    def _on_sssp_installed(self, change):
        self._model.installing_sssp = not change["new"]
        self.sssp_installation.layout.display = "none" if change["new"] else "block"

    def _on_code_activation_change(self, change):
        self._toggle_code(change["owner"])

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
        self.code_widgets[code.name] = code_widget
        self._model.set_selected_codes()
        code.is_rendered = True

    def _on_code_selection_change(self, _):
        self._update_submission_blockers()

    def _submit(self):
        """Submit the work chain with the current inputs."""

        parameters = deepcopy(self._model.input_parameters)
        builder = self._create_builder(parameters)

        with self.hold_trait_notifications():
            process = aiida_submit(builder)

            process.label = self._model.process_label
            process.description = self._model.process_description
            # since AiiDA data node may exist in the ui_parameters,
            # we serialize it to yaml
            process.base.extras.set("ui_parameters", serialize(parameters))
            # store the workchain name in extras, this will help to filter the workchain in the future
            process.base.extras.set("workchain", parameters["workchain"])  # type: ignore
            process.base.extras.set(
                "structure",
                self._model.input_structure.get_formula(),
            )
            self._model.process = process

    def _update(self):
        with self.hold_trait_notifications():
            self._update_active_codes()
            self._update_process_label()
            self._update_submission_blockers()

    def _update_submission_blockers(self):
        self.internal_submission_blockers = list(self._check_submission_blockers())

    def _update_submission_blocker_message(self):
        blockers = self.internal_submission_blockers + self.external_submission_blockers
        if any(blockers):
            fmt_list = "\n".join(f"<li>{item}</li>" for item in sorted(blockers))
            self._model.submission_blocker_messages = f"""
                <div class="alert alert-info">
                    <b>The submission is blocked due to the following reason(s):</b>
                    <ul>
                        {fmt_list}
                    </ul>
                </div>
            """
        else:
            self._model.submission_blocker_messages = ""

    def _check_submission_blockers(self):
        """Validate the resource inputs and identify blockers for the submission."""

        # Do not submit while any of the background setup processes are running.
        if self._model.installing_qe or self._model.installing_sssp:
            yield "Background setup processes must finish."

        # SSSP library not installed
        if not self._model.sssp_installed:
            yield "The SSSP library is not installed."

        # No pw code selected (this is ignored while the setup process is running).
        pw_code = self._model.get_code(identifier="dft", name="pw")
        if pw_code and not pw_code.selected and not self._model.installing_qe:
            yield ("No pw code selected")

        # code related to the selected property is not installed
        properties = self._model.get_properties()
        message = "Calculating the {property} property requires code {code} to be set."
        for identifier, codes in self._model.get_codes():
            if identifier in properties:
                for code in codes.values():
                    if not code.is_ready:
                        yield message.format(property=identifier, code=code.description)

        # check if the QEAppComputationalResourcesWidget is used
        for name, code in self._model.get_codes(flat=True):
            # skip if the code is not displayed, convenient for the plugin developer
            if not code.is_ready:
                continue
            if not issubclass(
                code.setup_widget_class, QEAppComputationalResourcesWidget
            ):
                yield (
                    f"Error: hi, plugin developer, please use the QEAppComputationalResourcesWidget from aiidalab_qe.common.widgets for code {name}."
                )

    def _update_active_codes(self):
        """Hide code if no related property is selected."""
        for name, code in self._model.get_codes(flat=True):
            if name != "pw":
                code.deactivate()
        properties = self._model.get_properties()
        for identifier, codes in self._model.get_codes():
            if identifier in properties:
                for code in codes.values():
                    code.activate()

    def _update_process_label(self):
        """Generate a label for the work chain based on the input parameters."""
        if not self._model.input_structure:
            return
        formula = self._model.input_structure.get_formula()
        workchain_data = self._model.input_parameters.get(
            "workchain",
            {"properties": []},
        )
        properties = [p for p in workchain_data["properties"] if p != "relax"]
        relax_type = workchain_data.get("relax_type", "none")
        if relax_type != "none":
            relax_info = "structure is relaxed"
        else:
            relax_info = "structure is not relaxed"
        if not properties:
            properties_info = ""
        else:
            properties_info = f", properties on {', '.join(properties)}"

        label = f"{formula} {relax_info} {properties_info}"
        self._model.process_label = label

    def _create_builder(self, parameters) -> ProcessBuilderNamespace:
        """Create the builder for the `QeAppWorkChain` submit."""

        submission_parameters = self._get_submission_parameters()
        parameters.update(submission_parameters)

        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self._model.input_structure,
            parameters=deepcopy(parameters),  # TODO why deepcopy again?
        )

        codes = submission_parameters["codes"]

        builder.relax.base.pw.metadata.options.resources = {
            "num_machines": codes.get("pw")["nodes"],
            "num_mpiprocs_per_machine": codes.get("pw")["ntasks_per_node"],
            "num_cores_per_mpiproc": codes.get("pw")["cpus_per_task"],
        }
        builder.relax.base.pw.metadata.options["max_wallclock_seconds"] = codes.get(
            "pw"
        )["max_wallclock_seconds"]
        builder.relax.base.pw.parallelization = orm.Dict(
            dict=codes["pw"]["parallelization"]
        )

        return builder

    def _get_submission_parameters(self):
        submission_parameters = self._model.get_model_state()
        for name, code_widget in self.code_widgets.items():
            if name in submission_parameters["codes"]:
                for key, value in code_widget.parameters.items():
                    if key != "code":
                        submission_parameters["codes"][name][key] = value
        return submission_parameters

    def _update_state(self, _=None):
        if self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
        elif self.previous_step_state is not self.State.SUCCESS:
            self.state = self.State.INIT
        elif self._model.process is not None:
            self.state = self.State.SUCCESS
        else:
            if self.is_blocked:
                self.state = self.State.READY
            else:
                self.internal_submission_blockers = []
                self.state = self.state.CONFIGURED
