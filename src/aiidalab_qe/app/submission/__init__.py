# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""
from __future__ import annotations

import os
import typing as t

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiida.common import NotExistent
from aiida.engine import ProcessBuilderNamespace, submit
from aiidalab_widgets_base import ComputationalResourcesWidget, WizardAppWidgetStep
from IPython.display import display

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.workflows import QeAppWorkChain

from .resource import ParallelizationSettings, ResourceSelectionWidget
from .sssp import SSSPInstallWidget

PROTOCOL_PSEUDO_MAP = {
    "fast": "SSSP/1.2/PBE/efficiency",
    "moderate": "SSSP/1.2/PBE/efficiency",
    "precise": "SSSP/1.2/PBE/precision",
}


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
    _submission_blockers = tl.List(tl.Unicode())

    def __init__(self, qe_auto_setup=True, **kwargs):
        self.message_area = ipw.Output()
        self._submission_blocker_messages = ipw.HTML()

        self.pw_code = ComputationalResourcesWidget(
            description="pw.x:", default_calc_job_plugin="quantumespresso.pw"
        )
        self.dos_code = ComputationalResourcesWidget(
            description="dos.x:",
            default_calc_job_plugin="quantumespresso.dos",
        )
        self.projwfc_code = ComputationalResourcesWidget(
            description="projwfc.x:",
            default_calc_job_plugin="quantumespresso.projwfc",
        )

        self.resources_config = ResourceSelectionWidget()
        self.parallelization = ParallelizationSettings()

        self.set_selected_codes(DEFAULT_PARAMETERS)
        self.set_resource_defaults()

        self.pw_code.observe(self._update_state, "value")
        self.pw_code.observe(self._update_resources, "value")
        self.dos_code.observe(self._update_state, "value")
        self.projwfc_code.observe(self._update_state, "value")

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
        self.sssp_installation_status = SSSPInstallWidget(auto_start=qe_auto_setup)
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

        super().__init__(
            children=[
                self.codes_title,
                self.codes_help,
                self.pw_code,
                self.dos_code,
                self.projwfc_code,
                self.resources_config,
                self.parallelization,
                self.message_area,
                self.sssp_installation_status,
                self.qe_setup_status,
                self._submission_blocker_messages,
                self.submit_button,
            ]
        )

    @tl.observe("_submission_blockers")
    def _observe_submission_blockers(self, change):
        if change["new"]:
            fmt_list = "\n".join((f"<li>{item}</li>" for item in sorted(change["new"])))
            self._submission_blocker_messages.value = f"""
                <div class="alert alert-info">
                <strong>The submission is blocked, due to the following reason(s):</strong>
                <ul>{fmt_list}</ul></div>"""
        else:
            self._submission_blocker_messages.value = ""

    def _identify_submission_blockers(self):
        # Do not submit while any of the background setup processes are running.
        if self.qe_setup_status.busy or self.sssp_installation_status.busy:
            yield "Background setup processes must finish."

        # No code selected (this is ignored while the setup process is running).
        if self.pw_code.value is None and not self.qe_setup_status.busy:
            yield ("No pw code selected")

        # No code selected for pdos (this is ignored while the setup process is running).
        if (
            "pdos" in self.input_parameters.get("workchain", {}).get("properties", [])
            and (self.dos_code.value is None or self.projwfc_code.value is None)
            and not self.qe_setup_status.busy
        ):
            yield "Calculating the PDOS requires both dos.x and projwfc.x to be set."

        # SSSP library not installed
        if not self.sssp_installation_status.installed:
            yield "The SSSP library is not installed."

        if (
            "pdos" in self.input_parameters.get("workchain", {}).get("properties", [])
            and not any(
                [
                    self.pw_code.value is None,
                    self.dos_code.value is None,
                    self.projwfc_code.value is None,
                ]
            )
            and len(
                set(
                    (
                        orm.load_code(self.pw_code.value).computer.pk,
                        orm.load_code(self.dos_code.value).computer.pk,
                        orm.load_code(self.projwfc_code.value).computer.pk,
                    )
                )
            )
            != 1
        ):
            yield (
                "All selected codes must be installed on the same computer. This is because the "
                "PDOS calculations rely on large files that are not retrieved by AiiDA."
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
            self._submission_blockers = blockers
            self.state = self.State.READY
            return

        self._submission_blockers = []
        self.state = self.state.CONFIGURED

    def _toggle_install_widgets(self, change):
        if change["new"]:
            self.children = [
                child for child in self.children if child is not change["owner"]
            ]

    def _auto_select_code(self, change):
        if change["new"] and not change["old"]:
            for code in [
                "pw_code",
                "dos_code",
                "projwfc_code",
            ]:
                try:
                    code_widget = getattr(self, code)
                    code_widget.refresh()
                    code_widget.value = orm.load_code(DEFAULT_PARAMETERS[code]).uuid
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

    def _update_resources(self, change):
        if change["new"] and (
            change["old"] is None
            or orm.load_code(change["new"]).computer.pk
            != orm.load_code(change["old"]).computer.pk
        ):
            self.set_resource_defaults(orm.load_code(change["new"]).computer)

    def set_resource_defaults(self, computer=None):
        if computer is None or computer.hostname == "localhost":
            self.resources_config.num_nodes.disabled = True
            self.resources_config.num_nodes.value = 1
            self.resources_config.num_cpus.max = os.cpu_count()
            self.resources_config.num_cpus.value = 1
            self.resources_config.num_cpus.description = "CPUs"
            self.parallelization.npools.value = 1
        else:
            default_mpiprocs = computer.get_default_mpiprocs_per_machine()
            self.resources_config.num_nodes.disabled = False
            self.resources_config.num_cpus.max = default_mpiprocs
            self.resources_config.num_cpus.value = default_mpiprocs
            self.resources_config.num_cpus.description = "CPUs/node"
            self.parallelization.npools.value = self._get_default_parallelization()

        self._check_resources()

    def _get_default_parallelization(self):
        """A _very_ rudimentary approach for obtaining a minimal npools setting."""
        num_mpiprocs = (
            self.resources_config.num_nodes.value * self.resources_config.num_cpus.value
        )

        for i in range(1, num_mpiprocs + 1):
            if num_mpiprocs % i == 0 and num_mpiprocs // i < self.MAX_MPI_PER_POOL:
                return i

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

    @tl.observe("previous_step_state")
    def _observe_input_structure(self, _):
        self._update_state()
        self.set_pdos_status()

    @tl.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                self.input_structure = process_node.inputs.structure
                builder_parameters = process_node.base.extras.get(
                    "builder_parameters", None
                )
                if builder_parameters is not None:
                    self.set_selected_codes(builder_parameters)
            self._update_state()

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.submit()

    def set_selected_codes(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        # Codes
        def _get_code_uuid(code):
            if code is not None:
                try:
                    return orm.load_code(code).uuid
                except NotExistent:
                    return None

        with self.hold_trait_notifications():
            # Codes
            self.pw_code.value = _get_code_uuid(parameters["pw_code"])
            self.dos_code.value = _get_code_uuid(parameters["dos_code"])
            self.projwfc_code.value = _get_code_uuid(parameters["projwfc_code"])

    def set_pdos_status(self):
        if "pdos" in self.input_parameters.get("workchain", {}).get("properties", []):
            self.dos_code.code_select_dropdown.disabled = False
            self.projwfc_code.code_select_dropdown.disabled = False
        else:
            self.dos_code.code_select_dropdown.disabled = True
            self.projwfc_code.code_select_dropdown.disabled = True

    def submit(self, _=None):
        """Submit the work chain with the current inputs."""
        builder = self._create_builder()
        extra_parameters = self._create_extra_report_parameters()

        with self.hold_trait_notifications():
            process = submit(builder)

            # Set the builder parameters on the work chain
            builder_parameters = self._extract_report_parameters(
                builder, extra_parameters
            )
            process.label = self._generate_label()
            process.base.extras.set("builder_parameters", builder_parameters)
            self.process = process

        self._update_state()

    def _generate_label(self) -> dict:
        """Generate a label for the work chain based on the input parameters."""
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
            properties_info = f"properties on {', '.join(properties)}"

        label = "{} {} {}".format(formula, relax_info, properties_info)
        return label

    def _create_builder(self) -> ProcessBuilderNamespace:
        """Create the builder for the `QeAppWorkChain` submit."""
        from copy import deepcopy

        pw_code = self.pw_code.value
        dos_code = self.dos_code.value
        projwfc_code = self.projwfc_code.value

        parameters = deepcopy(self.input_parameters)
        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            pw_code=orm.load_code(pw_code),
            dos_code=orm.load_code(dos_code),
            projwfc_code=orm.load_code(projwfc_code),
            parameters=parameters,
        )

        resources = {
            "num_machines": self.resources_config.num_nodes.value,
            "num_mpiprocs_per_machine": self.resources_config.num_cpus.value,
        }

        npool = self.parallelization.npools.value
        self._update_builder(builder, resources, npool, self.MAX_MPI_PER_POOL)

        return builder

    def _update_builder(self, buildy, resources, npools, max_mpi_per_pool):
        """Update the resources and parallelization of the ``QeAppWorkChain`` builder."""
        for k, v in buildy.items():
            if isinstance(v, (dict, ProcessBuilderNamespace)):
                if k == "pw" and v["pseudos"]:
                    v["parallelization"] = orm.Dict(dict={"npool": npools})
                if k == "projwfc":
                    v["settings"] = orm.Dict(dict={"cmdline": ["-nk", str(npools)]})
                if k == "dos":
                    v["metadata"]["options"]["resources"] = {
                        "num_machines": 1,
                        "num_mpiprocs_per_machine": min(
                            max_mpi_per_pool,
                            resources["num_mpiprocs_per_machine"],
                        ),
                    }
                    # Continue to the next item to avoid overriding the resources in the
                    # recursive `update_builder` call.
                    continue
                if k == "resources":
                    buildy["resources"] = resources
                else:
                    self._update_builder(v, resources, npools, max_mpi_per_pool)

    def _create_extra_report_parameters(self) -> dict[str, t.Any]:
        """This method will also create a dictionary of the parameters that were not
        readably represented in the builder, which will be used to the report.
        It is stored in the `extra_report_parameters`.
        """
        input_parameters = self.input_parameters

        # Construct the extra report parameters needed for the report
        extra_report_parameters = {
            "relax_type": input_parameters["workchain"]["relax_type"],
            "electronic_type": input_parameters["workchain"]["electronic_type"],
            "spin_type": input_parameters["workchain"]["spin_type"],
            "protocol": input_parameters["workchain"]["protocol"],
            "initial_magnetic_moments": input_parameters["advanced"][
                "initial_magnetic_moments"
            ],
            "properties": input_parameters["workchain"]["properties"],
        }

        # Periodicity
        periodicity_mapping = {
            (True, True, True): "xyz",
            (True, True, False): "xy",
            (True, False, False): "x",
        }

        periodicity = periodicity_mapping.get(self.input_structure.pbc, "xyz")
        extra_report_parameters.update({"periodicity": periodicity})

        # update pseudo family information to extra_report_parameters
        pseudo_family = input_parameters["advanced"]["pseudo_family"]
        pseudo_family = PROTOCOL_PSEUDO_MAP[input_parameters["workchain"]["protocol"]]

        pseudo_family_info = pseudo_family.split("/")
        if pseudo_family_info[0] == "SSSP":
            pseudo_protocol = pseudo_family_info[3]
        elif pseudo_family_info[0] == "PseudoDojo":
            pseudo_protocol = pseudo_family_info[4]
        extra_report_parameters.update(
            {
                "pseudo_family": pseudo_family,
                "pseudo_library": pseudo_family_info[0],
                "pseudo_version": pseudo_family_info[1],
                "functional": pseudo_family_info[2],
                "pseudo_protocol": pseudo_protocol,
            }
        )

        # store codes info into extra_report_parameters for loading the process
        pw_code = self.pw_code.value
        dos_code = self.dos_code.value
        projwfc_code = self.projwfc_code.value

        extra_report_parameters.update(
            {
                "pw_code": pw_code,
                "dos_code": dos_code,
                "projwfc_code": projwfc_code,
            }
        )

        return extra_report_parameters

    @staticmethod
    def _extract_report_parameters(
        builder, extra_report_parameters
    ) -> dict[str, t.Any]:
        """Extract (recover) the parameters for report from the builder.

        There are some parameters that are not stored in the builder, but can be extracted
        directly from the widgets, such as the ``pseudo_family`` and ``relax_type``.
        """
        parameters = {
            "run_relax": "relax" in builder.properties,
            "run_bands": "bands" in builder.properties,
            "run_pdos": "pdos" in builder.properties,
        }

        # Extract the pw calculation parameters from the builder

        # energy_cutoff is same for all pw calculations when pseudopotentials are fixed
        # as well as the smearing settings (semaring and degauss) and scf kpoints distance
        # read from the first pw calculation of relax workflow.
        # It is safe then to extract these parameters from the first pw calculation, since the
        # builder is anyway set with subworkchain inputs even it is not run which controlled by
        # the properties inputs.
        energy_cutoff_wfc = builder.relax.base["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
        energy_cutoff_rho = builder.relax.base["pw"]["parameters"]["SYSTEM"]["ecutrho"]
        occupation = builder.relax.base["pw"]["parameters"]["SYSTEM"]["occupations"]
        scf_kpoints_distance = builder.relax.base.kpoints_distance.value

        parameters.update(
            {
                "energy_cutoff_wfc": energy_cutoff_wfc,
                "energy_cutoff_rho": energy_cutoff_rho,
                "occupation": occupation,
                "scf_kpoints_distance": scf_kpoints_distance,
            }
        )

        if occupation == "smearing":
            parameters["degauss"] = builder.relax.base["pw"]["parameters"]["SYSTEM"][
                "degauss"
            ]
            parameters["smearing"] = builder.relax.base["pw"]["parameters"]["SYSTEM"][
                "smearing"
            ]

        if "bands" in builder.properties:
            parameters[
                "bands_kpoints_distance"
            ] = builder.bands.bands_kpoints_distance.value

        parameters["tot_charge"] = builder.relax.base["pw"]["parameters"]["SYSTEM"].get(
            "tot_charge", 0.0
        )

        # parameters from extra_report_parameters
        for k, v in extra_report_parameters.items():
            parameters.update({k: v})

        return parameters

    def reset(self):
        with self.hold_trait_notifications():
            self.process = None
            self.input_structure = None
