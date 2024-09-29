"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl
from IPython.display import display

from aiida import orm
from aiida.engine import ProcessBuilderNamespace
from aiida.engine import submit as aiida_submit
from aiida.orm.utils.serialize import serialize
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.panel import CodesDict
from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.common.setup_pseudos import PseudosInstallWidget
from aiidalab_qe.common.widgets import (
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)
from aiidalab_qe.workflows import QeAppWorkChain
from aiidalab_widgets_base import WizardAppWidgetStep

from .model import SubmissionModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    process = tl.Instance(orm.WorkChainNode, allow_none=True)
    input_parameters = tl.Dict()

    internal_submission_blockers = tl.List(tl.Unicode())
    external_submission_blockers = tl.List(tl.Unicode())

    def __init__(self, model: SubmissionModel, qe_auto_setup=True, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workflow submission panel")],
            **kwargs,
        )

        self._model = model

        # # TODO for testing only - remove in PR
        # self._model.observe(
        #     lambda change: print(change["new"]),
        #     "input_parameters",
        # )

        self.qe_auto_setup = qe_auto_setup

        # TODO consider restructuring model's codes trait to remove this
        self.code_fetchers: dict[str, CodesDict] = get_entry_items(
            "aiidalab_qe.properties", "code"
        )

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

        self.code_children = []

        # add default pw code
        pw_code = PwCodeResourceSetupWidget(
            description="pw.x:",
            default_calc_job_plugin="quantumespresso.pw",
        )
        pw_code.observe(
            self._update_state,
            "value",
        )
        self.code_children.append(pw_code)
        self._model.add_code("pw", pw_code)

        for codes_fetcher in self.code_fetchers.values():
            for name, code in codes_fetcher.items():
                self._model.add_code(name, code)
                code.observe(
                    self._update_state,
                    "value",
                )
                self.code_children.append(code)

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

        self.sssp_installation = PseudosInstallWidget(auto_start=self.qe_auto_setup)
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

        self.qe_setup = QESetupWidget(auto_start=self.qe_auto_setup)
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
            *self.code_children,
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

        # set default codes
        self._model.set_selected_codes(DEFAULT["codes"])

        with self.hold_trait_notifications():
            ipw.dlink(
                (self._model, "input_structure"),
                (self, "input_structure"),
            )
            ipw.dlink(
                (self._model, "process"),
                (self, "process"),
            )
            ipw.dlink(
                (self._model, "input_parameters"),
                (self, "input_parameters"),
            )

        self.rendered = True

    # TODO seems to only be setting codes? rename for clarity? same with setter
    def get_submission_parameters(self):
        return {
            "codes": self._model.get_selected_codes(),
        }

    def set_submission_parameters(self, parameters):
        self._model.set_model_state(parameters)

    def set_process(self, process):
        self._model.process = process

    def submit(self, _=None):
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
            process.base.extras.set("structure", self.input_structure.get_formula())
            self._model.process = process

        self._update_state()

    def reset(self):
        with self.hold_trait_notifications():
            self._model.reset()
            self._model.set_selected_codes(DEFAULT["codes"])

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update()

    @tl.observe("input_structure")
    def _on_input_structure_change(self, _):
        self._update()

    @tl.observe("process")
    def _on_process_change(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            # TODO why here? Do we not populate traits earlier that would cover this?
            if process_node is not None:
                self._model.input_structure = process_node.inputs.structure
            self._update_state()

    @tl.observe("input_parameters")
    def _on_input_parameters_change(self, _):
        self._update()

    @tl.observe("internal_submission_blockers", "external_submission_blockers")
    def _on_submission_blockers_change(self, _):
        self._update_submission_blocker_message()

    def _on_submission(self, _):
        self.submit_button.disabled = True
        self.submit()

    def _on_installation_change(self, _):
        with self.hold_trait_notifications():
            self._update_submission_blockers()
            self._update_state()

    def _on_qe_installed(self, change):
        self._model.installing_qe = not change["new"]
        self.qe_setup.layout.display = "none" if change["new"] else "block"
        self._auto_select_code(change)

    def _on_sssp_installed(self, change):
        self._model.installing_sssp = not change["new"]
        self.sssp_installation.layout.display = "none" if change["new"] else "block"

    def _update(self):
        with self.hold_trait_notifications():
            self._update_submission_blockers()
            self._update_codes_display()
            self._update_process_label()
            self._update_state()

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

        # No pw code selected (this is ignored while the setup process is running).
        pw_code = self._model.get_code("pw")
        if pw_code and pw_code.value is None and not self._model.installing_qe:
            yield ("No pw code selected")

        # code related to the selected property is not installed
        properties = self._model.input_parameters.get("workchain", {}).get(
            "properties", []
        )
        for identifier in properties:
            codes_fetcher = self.code_fetchers.get(identifier, {})
            for name, code in codes_fetcher.items():
                if code.value is None:
                    yield f"Calculating the {identifier} property requires code {name} to be set."

        # SSSP library not installed
        if not self._model.sssp_installed:
            yield "The SSSP library is not installed."

        # check if the QEAppComputationalResourcesWidget is used
        for name, code in self._model.codes.items():
            # skip if the code is not displayed, convenient for the plugin developer
            if code.layout.display == "none":
                continue
            if not isinstance(code, QEAppComputationalResourcesWidget):
                yield (
                    f"Error: hi, plugin developer, please use the QEAppComputationalResourcesWidget from aiidalab_qe.common.widgets for code {name}."
                )

    def _auto_select_code(self, change):
        if change["new"] and not change["old"]:
            self._model.set_selected_codes(DEFAULT["codes"])

    def _show_alert_message(self, message, alert_class="info"):
        with self.message_area:
            display(
                ipw.HTML(
                    self._ALERT_MESSAGE.format(
                        alert_class=alert_class,
                        message=message,
                    )
                )
            )

    def _update_codes_display(self):
        """Hide code if no related property is selected."""
        for name, code in self._model.codes.items():
            if name == "pw":
                continue
            code.layout.display = "none"
        for _, code in self._model.codes.items():
            code.layout.display = "block"

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

        submission_parameters = self.get_submission_parameters()
        parameters.update(submission_parameters)
        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
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
