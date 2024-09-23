"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl
from IPython.display import display

from aiida import orm
from aiida.common import NotExistent
from aiida.engine import ProcessBuilderNamespace
from aiida.engine import submit as aiida_submit
from aiida.orm.utils.serialize import serialize
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
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

from .model import SubmissionModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    # This number provides a rough estimate for how many MPI tasks are needed
    # for a given structure.
    NUM_SITES_PER_MPI_TASK_DEFAULT = 6

    # Warn the user if they are trying to run calculations for a large
    # structure on localhost.
    RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD = 10

    # Put a limit on how many MPI tasks you want to run per k-pool by default
    MAX_MPI_PER_POOL = 20

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

        # TODO for testing only - remove in PR
        self._model.observe(
            lambda change: print(change["new"]),
            "input_parameters",
        )

        self.qe_auto_setup = qe_auto_setup
        self._submission_blocker_messages = ipw.HTML()

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

        self.message_area = ipw.Output()

        self.pw_code = PwCodeResourceSetupWidget(
            description="pw.x:",
            default_calc_job_plugin="quantumespresso.pw",
        )

        self.pw_code.observe(
            self._update_state,
            "value",
        )

        self.codes = {"pw": self.pw_code}
        self.code_children = [
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
            self.pw_code,
        ]

        # add plugin's entry points
        self.code_entries = get_entry_items("aiidalab_qe.properties", "code")
        for _, entry_point in self.code_entries.items():
            for identifier, code in entry_point().items():
                self.codes[identifier] = code
                code.observe(
                    self._update_state,
                    "value",
                )
                self.code_children.append(self.codes[identifier])

        # set process label and description
        self.process_label = ipw.Text(
            description="Label:",
            layout=ipw.Layout(width="auto", indent="0px"),
        )
        self.process_description = ipw.Textarea(
            description="Description",
            layout=ipw.Layout(width="auto", indent="0px"),
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

        self.submit_button.on_click(self._on_submit_button_clicked)

        self.sssp_installation_status = PseudosInstallWidget(
            auto_start=self.qe_auto_setup
        )
        self.sssp_installation_status.observe(
            self._update_state,
            ["busy", "installed"],
        )
        self.sssp_installation_status.observe(
            self._toggle_install_widgets,
            "installed",
        )

        self.qe_setup_status = QESetupWidget(auto_start=self.qe_auto_setup)
        self.qe_setup_status.observe(
            self._update_state,
            "busy",
        )
        self.qe_setup_status.observe(
            self._toggle_install_widgets,
            "installed",
        )
        self.qe_setup_status.observe(
            self._auto_select_code,
            "installed",
        )

        self.ui_parameters = {}

        self.children = [
            *self.code_children,
            self.message_area,
            self.sssp_installation_status,
            self.qe_setup_status,
            self._submission_blocker_messages,
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
        self._set_selected_codes(DEFAULT["codes"])

        self.observe(self._on_input_structure_change, "input_parameters")

        with self.hold_trait_notifications():
            ipw.dlink(
                (self._model, "input_structure"),
                (self, "input_structure"),
            )
            ipw.dlink(
                (self._model, "input_parameters"),
                (self, "input_parameters"),
            )
            ipw.dlink(
                (self, "process"),
                (self._model, "process"),
            )

        self.rendered = True

    def get_submission_parameters(self):
        """Get the parameters for the submission step."""
        return {
            "codes": self._get_selected_codes(),
        }

    def set_submission_parameters(self, parameters):
        # backward compatibility for v2023.11
        # which have a separate "resources" section for pw code
        if "resources" in parameters:
            parameters["codes"] = {
                key: {"code": value} for key, value in parameters["codes"].items()
            }
            parameters["codes"]["pw"]["nodes"] = parameters["resources"]["num_machines"]
            parameters["codes"]["pw"]["cpus"] = parameters["resources"][
                "num_mpiprocs_per_machine"
            ]
            parameters["codes"]["pw"]["parallelization"] = {
                "npool": parameters["resources"]["npools"]
            }
        self._set_selected_codes(parameters["codes"])
        # label and description are not stored in the parameters, but in the process directly
        if self.process:
            self.process_label.value = self.process.label
            self.process_description.value = self.process.description

    def submit(self, _=None):
        """Submit the work chain with the current inputs."""

        builder = self._create_builder()

        with self.hold_trait_notifications():
            process = aiida_submit(builder)

            process.label = self.process_label.value
            process.description = self.process_description.value
            # since AiiDA data node may exist in the ui_parameters,
            # we serialize it to yaml
            process.base.extras.set("ui_parameters", serialize(self.ui_parameters))
            # store the workchain name in extras, this will help to filter the workchain in the future
            process.base.extras.set("workchain", self.ui_parameters["workchain"])
            process.base.extras.set("structure", self.input_structure.get_formula())
            self.process = process

        self._update_state()

    def reset(self):
        """Reset the widget to its initial state."""
        with self.hold_trait_notifications():
            self.process = None
            self.input_structure = None
            self._set_selected_codes(DEFAULT["codes"])

    @tl.observe("input_structure")
    def _on_input_structure_change(self, _):
        self._update_state()
        self._update_codes_display()
        self._update_process_label()

    @tl.observe("process")
    def _on_process_change(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            # TODO why here? Do we not populate traits earlier that would cover this?
            if process_node is not None:
                self.input_structure = process_node.inputs.structure
            self._update_state()

    @tl.observe("internal_submission_blockers", "external_submission_blockers")
    def _on_submission_blockers_change(self, _):
        self._update_submission_blocker_message()

    def _update_submission_blocker_message(self):
        blockers = self.internal_submission_blockers + self.external_submission_blockers
        if any(blockers):
            fmt_list = "\n".join(f"<li>{item}</li>" for item in sorted(blockers))
            self._submission_blocker_messages.value = f"""
                <div class="alert alert-info">
                    <b>The submission is blocked, due to the following reason(s):</b>
                    <ul>
                        {fmt_list}
                    </ul>
                </div>
            """
        else:
            self._submission_blocker_messages.value = ""

    def _update_state(self, _=None):
        # If the previous step has failed, this should fail as well.
        if self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
            return
        # Do not interact with the user if they haven't successfully completed the previous step.
        elif self.previous_step_state is not self.State.SUCCESS:
            self.state = self.State.INIT
            return

        # Process is already running.
        if self.process is not None:
            self.state = self.State.SUCCESS
            return

        blockers = list(self._identify_submission_blockers())
        if any(blockers):
            self.internal_submission_blockers = blockers
            self.state = self.State.READY
            return

        self.internal_submission_blockers = []
        self.state = self.state.CONFIGURED

    def _identify_submission_blockers(self):
        """Validate the resource inputs and identify blockers for the submission."""

        # Do not submit while any of the background setup processes are running.
        if self.qe_setup_status.busy or self.sssp_installation_status.busy:
            yield "Background setup processes must finish."

        # No pw code selected (this is ignored while the setup process is running).
        if self.pw_code.value is None and not self.qe_setup_status.busy:
            yield ("No pw code selected")
        # code related to the selected property is not installed
        properties = self.input_parameters.get("workchain", {}).get("properties", [])
        for identifer in properties:
            for name, code in self.code_entries.get(identifer, {}).items():
                if code.value is None:
                    yield f"Calculating the {identifer} property requires code {name} to be set."
        # SSSP library not installed
        if not self.sssp_installation_status.installed:
            yield "The SSSP library is not installed."

        # check if the QEAppComputationalResourcesWidget is used
        for name, code in self.codes.items():
            # skip if the code is not displayed, convenient for the plugin developer
            if code.layout.display == "none":
                continue
            if not isinstance(code, QEAppComputationalResourcesWidget):
                yield (
                    f"Error: hi, plugin developer, please use the QEAppComputationalResourcesWidget from aiidalab_qe.common.widgets for code {name}."
                )

    def _get_selected_codes(self):
        """Get the codes selected in the GUI.

        return: A dict with the code names as keys and the code UUIDs as values.
        """
        codes = {
            key: code.parameters
            for key, code in self.codes.items()
            if code.layout.display != "none"
        }
        return codes

    def _set_selected_codes(self, code_data):
        """Set the inputs in the GUI based on a set of codes."""

        # Codes
        def _get_code_uuid(code):
            if code is not None:
                try:
                    return orm.load_code(code).uuid
                except NotExistent:
                    return None

        with self.hold_trait_notifications():
            for name, code in self.codes.items():
                if name not in code_data:
                    continue
                # check if the code is installed and usable
                # note: if code is imported from another user, it is not usable and thus will not be
                # treated as an option in the ComputationalResourcesWidget.
                code_options = [
                    o[1] for o in code.code_selection.code_select_dropdown.options
                ]
                if _get_code_uuid(code_data.get(name)["code"]) in code_options:
                    # get code uuid from code label in case of using DEFAULT_PARAMETERS
                    code_data.get(name)["code"] = _get_code_uuid(
                        code_data.get(name)["code"]
                    )
                    code.parameters = code_data.get(name)

    def _toggle_install_widgets(self, change):
        if change["new"]:
            self.children = [
                child for child in self.children if child is not change["owner"]
            ]

    def _auto_select_code(self, change):
        if change["new"] and not change["old"]:
            self._set_selected_codes(DEFAULT["codes"])

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

    def _check_resources(self):
        """Check whether the currently selected resources will be sufficient and warn if not."""
        if not self.pw_code.value:
            return  # No code selected, nothing to do.

        num_cpus = self.resources_config.num_cpus.value
        on_localhost = (
            orm.load_node(self.pw_code.value).computer.hostname == "localhost"
        )
        if self.pw_code.value and on_localhost and num_cpus > 1:
            self._show_alert_message(
                "The selected code would be executed on the local host, but "
                "the number of CPUs is larger than one. Please review "
                "the configuration and consider to select a code that runs "
                "on a larger system if necessary.",
                alert_class="warning",
            )
        elif (
            self.input_structure
            and on_localhost
            and len(self.input_structure.sites)
            > self.RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD
        ):
            self._show_alert_message(
                "The selected code would be executed on the local host, but the "
                "number of sites of the selected structure is relatively large. "
                "Consider to select a code that runs on a larger system if "
                "necessary.",
                alert_class="warning",
            )

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.submit()

    def _update_codes_display(self):
        """Hide code if no related property is selected."""
        # hide all codes except pw
        for name, code in self.codes.items():
            if name == "pw":
                continue
            code.layout.display = "none"
        properties = self.input_parameters.get("workchain", {}).get("properties", [])
        # show the code if the related property is selected.
        for identifier in properties:
            for code in self.code_entries.get(identifier, {}).values():
                code.layout.display = "block"

    def _update_process_label(self) -> dict:
        """Generate a label for the work chain based on the input parameters."""
        if not self.input_structure:
            return ""
        formula = self.input_structure.get_formula()
        workchain_data = self.input_parameters.get("workchain", {"properties": []})
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
        self.process_label.value = label

    def _create_builder(self) -> ProcessBuilderNamespace:
        """Create the builder for the `QeAppWorkChain` submit."""

        self.ui_parameters = deepcopy(self.input_parameters)
        # add codes and resource info into ui_parameters
        submission_parameters = self.get_submission_parameters()
        self.ui_parameters.update(submission_parameters)
        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            parameters=deepcopy(self.ui_parameters),
        )

        self._update_builder(builder, submission_parameters["codes"])

        return builder

    def _update_builder(self, builder, codes):
        """Update the resources and parallelization of the ``relax`` builder."""
        # update resources
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
