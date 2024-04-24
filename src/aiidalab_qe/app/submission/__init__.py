# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiida.common import NotExistent
from aiida.engine import ProcessBuilderNamespace, submit
from aiidalab_widgets_base import WizardAppWidgetStep
from IPython.display import display

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.common.setup_pseudos import PseudosInstallWidget
from aiidalab_qe.common.widgets import (
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)
from aiidalab_qe.workflows import QeAppWorkChain


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    codes_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Codes</h4></div>"""
    )
    codes_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom:
        10px"> Select the code to use for running the calculations. The codes
        on the local machine (localhost) are installed by default, but you can
        configure new ones on potentially more powerful machines by clicking on
        "Setup new code".</div>"""
    )
    process_label_help = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Labeling Your Job</h4>
        <p style="line-height: 140%; padding-top: 0px; padding-bottom:
        10px"> Label your job and provide a brief description. These details help identify the job later and make the search process easier. While optional, adding a description is recommended for better clarity.</p>
        </div>"""
    )

    # This number provides a rough estimate for how many MPI tasks are needed
    # for a given structure.
    NUM_SITES_PER_MPI_TASK_DEFAULT = 6

    # Warn the user if they are trying to run calculations for a large
    # structure on localhost.
    RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD = 10

    # Put a limit on how many MPI tasks you want to run per k-pool by default
    MAX_MPI_PER_POOL = 20

    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    process = tl.Instance(orm.WorkChainNode, allow_none=True)
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)
    input_parameters = tl.Dict()
    internal_submission_blockers = tl.List(tl.Unicode())
    external_submission_blockers = tl.List(tl.Unicode())

    def __init__(self, qe_auto_setup=True, **kwargs):
        self.message_area = ipw.Output()
        self._submission_blocker_messages = ipw.HTML()

        self.pw_code = PwCodeResourceSetupWidget(
            description="pw.x:", default_calc_job_plugin="quantumespresso.pw"
        )

        self.pw_code.observe(self._update_state, "value")

        # add plugin's entry points
        self.codes = {"pw": self.pw_code}
        self.code_children = [
            self.codes_title,
            self.codes_help,
            self.pw_code,
        ]
        self.code_entries = get_entry_items("aiidalab_qe.properties", "code")
        for _, entry_point in self.code_entries.items():
            for name, code in entry_point.items():
                self.codes[name] = code
                code.observe(self._update_state, "value")
                self.code_children.append(self.codes[name])
        # set process label and description
        self.process_label = ipw.Text(
            description="Label:", layout=ipw.Layout(width="auto", indent="0px")
        )
        self.process_description = ipw.Textarea(
            description="Description", layout=ipw.Layout(width="auto", indent="0px")
        )
        #
        self.submit_button = ipw.Button(
            description="Submit",
            tooltip="Submit the calculation with the selected parameters.",
            icon="play",
            button_style="success",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )

        self.submit_button.on_click(self._on_submit_button_clicked)

        # The SSSP installation status widget shows the installation status of
        # the SSSP pseudo potentials and triggers the installation in case that
        # they are not yet installed. The widget will remain in a "busy" state
        # in case that the installation was already triggered elsewhere, e.g.,
        # by the start up scripts.  The submission is blocked while the
        # potentials are not yet installed.
        self.sssp_installation_status = PseudosInstallWidget(auto_start=qe_auto_setup)
        self.sssp_installation_status.observe(self._update_state, ["busy", "installed"])
        self.sssp_installation_status.observe(self._toggle_install_widgets, "installed")

        # The QE setup widget checks whether there are codes that match specific
        # expected labels (e.g. "pw-7.2@localhost") and triggers both the
        # installation of QE into a dedicated conda environment and the setup of
        # the codes in case that they are not already configured.
        self.qe_setup_status = QESetupWidget(auto_start=qe_auto_setup)
        self.qe_setup_status.observe(self._update_state, "busy")
        self.qe_setup_status.observe(self._toggle_install_widgets, "installed")
        self.qe_setup_status.observe(self._auto_select_code, "installed")
        self.ui_parameters = {}

        super().__init__(
            children=[
                *self.code_children,
                self.message_area,
                self.sssp_installation_status,
                self.qe_setup_status,
                self._submission_blocker_messages,
                self.process_label_help,
                self.process_label,
                self.process_description,
                self.submit_button,
            ]
        )
        # set default codes
        self.set_selected_codes(DEFAULT_PARAMETERS["codes"])

    @tl.observe("internal_submission_blockers", "external_submission_blockers")
    def _observe_submission_blockers(self, _change):
        """Observe the submission blockers and update the message area."""
        blockers = self.internal_submission_blockers + self.external_submission_blockers
        if any(blockers):
            fmt_list = "\n".join((f"<li>{item}</li>" for item in sorted(blockers)))
            self._submission_blocker_messages.value = f"""
                <div class="alert alert-info">
                <strong>The submission is blocked, due to the following reason(s):</strong>
                <ul>{fmt_list}</ul></div>"""
        else:
            self._submission_blocker_messages.value = ""

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

    def _toggle_install_widgets(self, change):
        if change["new"]:
            self.children = [
                child for child in self.children if child is not change["owner"]
            ]

    def _auto_select_code(self, change):
        if change["new"] and not change["old"]:
            for name, code in self.codes.items():
                if not DEFAULT_PARAMETERS["codes"].get(name):
                    continue
                try:
                    code.code_selection.refresh()
                    code.value = orm.load_code(DEFAULT_PARAMETERS["codes"][name]).uuid
                except NotExistent:
                    pass

    _ALERT_MESSAGE = """
        <div class="alert alert-{alert_class} alert-dismissible">
        <a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>
        <span class="closebtn" onclick="this.parentElement.style.display='none';">&times;</span>
        <strong>{message}</strong>
        </div>"""

    def _show_alert_message(self, message, alert_class="info"):
        with self.message_area:
            display(
                ipw.HTML(
                    self._ALERT_MESSAGE.format(alert_class=alert_class, message=message)
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

    @tl.observe("state")
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.submit_button.disabled = change["new"] != self.State.CONFIGURED

    @tl.observe("previous_step_state", "input_parameters")
    def _observe_input_structure(self, _):
        self._update_state()
        self.update_codes_display()
        self._update_process_label()

    @tl.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                self.input_structure = process_node.inputs.structure
            self._update_state()

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.submit()

    def get_selected_codes(self):
        """Get the codes selected in the GUI.

        return: A dict with the code names as keys and the code UUIDs as values.
        """
        codes = {
            key: code.parameters
            for key, code in self.codes.items()
            if code.layout.display != "none"
        }
        return codes

    def set_selected_codes(self, code_data):
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

    def update_codes_display(self):
        """Hide code if no related property is selected."""
        # hide all codes except pw
        for name, code in self.codes.items():
            if name == "pw":
                continue
            code.layout.display = "none"
        properties = self.input_parameters.get("workchain", {}).get("properties", [])
        # show the code if the related property is selected.
        for identifer in properties:
            for code in self.code_entries.get(identifer, {}).values():
                code.layout.display = "block"

    def submit(self, _=None):
        """Submit the work chain with the current inputs."""
        from aiida.orm.utils.serialize import serialize

        builder = self._create_builder()

        with self.hold_trait_notifications():
            process = submit(builder)

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

    def _update_process_label(self) -> dict:
        """Generate a label for the work chain based on the input parameters."""
        if not self.input_structure:
            return ""
        formula = self.input_structure.get_formula()
        properties = [
            p for p in self.input_parameters["workchain"]["properties"] if p != "realx"
        ]
        relax_type = self.input_parameters["workchain"].get("relax_type")
        if relax_type != "none":
            relax_info = "structure is relaxed"
        else:
            relax_info = "structure is not relaxed"
        if not properties:
            properties_info = ""
        else:
            properties_info = f", properties on {', '.join(properties)}"

        label = "{} {} {}".format(formula, relax_info, properties_info)
        self.process_label.value = label

    def _create_builder(self) -> ProcessBuilderNamespace:
        """Create the builder for the `QeAppWorkChain` submit."""
        from copy import deepcopy

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
        builder.relax.base.pw.parallelization = orm.Dict(
            dict=codes["pw"]["parallelization"]
        )

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
        self.set_selected_codes(parameters["codes"])
        # label and description are not stored in the parameters, but in the process directly
        if self.process:
            self.process_label.value = self.process.label
            self.process_description.value = self.process.description

    def get_submission_parameters(self):
        """Get the parameters for the submission step."""
        return {
            "codes": self.get_selected_codes(),
        }

    def reset(self):
        """Reset the widget to its initial state."""
        with self.hold_trait_notifications():
            self.process = None
            self.input_structure = None
            self.set_selected_codes(DEFAULT_PARAMETERS["codes"])
