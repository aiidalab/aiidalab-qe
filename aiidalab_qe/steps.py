"""Widgets for the submission of bands work chains.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
from math import ceil

import ipywidgets as ipw
import traitlets
from aiida.common import NotExistent
from aiida.engine import ProcessState, submit
from aiida.orm import ProcessNode, WorkChainNode, load_code
from aiida.plugins import DataFactory
from aiida_quantumespresso.common.types import ElectronicType, RelaxType, SpinType
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_widgets_base import (
    CodeDropdown,
    ProcessMonitor,
    ProcessNodesTreeWidget,
    WizardAppWidgetStep,
)
from IPython.display import display

from aiidalab_qe.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.pseudos import PseudoFamilySelector
from aiidalab_qe.setup_codes import QESetupWidget
from aiidalab_qe.sssp import SSSPInstallWidget
from aiidalab_qe.widgets import NodeViewWidget, ResourceSelectionWidget
from aiidalab_qe_workchain import QeAppWorkChain

StructureData = DataFactory("structure")
Float = DataFactory("float")


def update_resources(builder, resources):
    for k, v in builder.items():
        if isinstance(v, dict):
            if k == "resources":
                builder["resources"].update(resources)
            else:
                update_resources(v, resources)


class WorkChainSettings(ipw.VBox):

    structure_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Structure</h4></div>"""
    )
    structure_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
        By default, the workflow will optimize the provided geometry. Select "Structure
        as is" if this is not desired. You can either optimize the atomic positions ("Full geometry")
        and unit cell, or atomic positions only ("Atomic positions"). </div>"""
    )
    materials_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Below you can indicate both if the material is magnetic and a metal or insulator. For now only ferromagnetic configurations are possible, since antiferromagnetism is more complicated to study automatically. If you're not sure whether your material is insulating, choose "Metal", since the corresponding settings usually also work quite well for insulators.
        </div>"""
    )

    properties_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Properties</h4></div>"""
    )
    properties_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
        The band structure workflow will
        automatically detect the default path in reciprocal space using the
        <a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">
        SeeK-path tool</a>.</div>"""
    )

    protocol_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Protocol</h4></div>"""
    )
    protocol_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
        The "moderate" protocol represents a balanced trade-off between
        accuracy and speed. Choose the "fast" protocol for a faster calculation
        with less precision and the "precise" protocol that provides more
        accuracy but will take longer.</div>"""
    )

    def __init__(self, **kwargs):

        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "none"),
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ],
            value="positions_cell",
        )

        # SpinType: magnetic properties of material
        self.spin_type = ipw.Dropdown(
            options=[("Non-magnetic", "none"), ("Ferromagnetic", "collinear")],
            value=DEFAULT_PARAMETERS["spin_type"],
            description="Magnetism:",
            style={"description_width": "initial"},
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.Dropdown(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            value=DEFAULT_PARAMETERS["electronic_type"],
            description="Electronic Type:",
            style={"description_width": "initial"},
        )

        # Checkbox to see if the band structure should be calculated
        self.bands_run = ipw.Checkbox(
            description="",
            tooltip="Calculate the electronic band structure.",
            indent=False,
            value=True,
            layout=ipw.Layout(max_width="10%"),
        )

        # Checkbox to see if the PDOS should be calculated
        self.pdos_run = ipw.Checkbox(
            description="",
            tooltip="Calculate the electronic PDOS.",
            indent=False,
            value=True,
            layout=ipw.Layout(max_width="10%"),
        )

        # Work chain protocol
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        super().__init__(
            children=[
                self.structure_title,
                self.structure_help,
                self.relax_type,
                self.materials_help,
                self.spin_type,
                self.electronic_type,
                self.properties_title,
                ipw.HTML("Select which properties to calculate:"),
                ipw.HBox(children=[ipw.HTML("<b>Band structure</b>"), self.bands_run]),
                ipw.HBox(
                    children=[
                        ipw.HTML("<b>Projected density of states</b>"),
                        self.pdos_run,
                    ]
                ),
                self.properties_help,
                self.protocol_title,
                ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
                self.workchain_protocol,
                self.protocol_help,
            ],
            **kwargs,
        )


class KpointSettings(ipw.VBox):

    kpoints_distance_description = ipw.HTML(
        """<p>
        The k-points mesh density is set by the chosen <b>protocol</b>.
        The value below represents the maximum distance between the k-points in each direction of reciprocal space.
        Untick the box to override the default.
    </p>"""
    )

    # The default of `kpoints_distance` must be linked to the `protocol`
    kpoints_distance_default = traitlets.Float(default_value=0.15)

    def __init__(self, **kwargs):

        self.override_protocol_kpoints = ipw.Checkbox(
            description="Override default k-points distance.",
            indent=False,
            value=False,
        )
        self.kpoints_distance = ipw.FloatText(
            value=self.kpoints_distance_default,
            step=0.05,
            description="K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override_protocol_kpoints, "value"),
            (self.kpoints_distance, "disabled"),
            lambda override: not override,
        )
        self.kpoints_distance.observe(self.set_kpoints_distance, "value")
        self.override_protocol_kpoints.observe(self.set_kpoints_distance, "value")
        self.observe(self.set_kpoints_distance, "kpoints_distance_default")

        super().__init__(
            children=[
                self.kpoints_distance_description,
                ipw.HBox([self.override_protocol_kpoints, self.kpoints_distance]),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    def set_kpoints_distance(self, _=None):
        self.kpoints_distance.value = (
            self.kpoints_distance.value
            if self.override_protocol_kpoints.value
            else self.kpoints_distance_default
        )


class CodeSettings(ipw.VBox):

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

    def __init__(self, **kwargs):

        self.pw = CodeDropdown(
            input_plugin="quantumespresso.pw",
            description="pw.x:",
            setup_code_params={
                "computer": "localhost",
                "description": "pw.x in AiiDAlab container.",
                "label": "pw-6.7",
                "input_plugin": "quantumespresso.pw",
                "remote_abs_path": "/usr/bin/pw.x",
            },
        )
        self.dos = CodeDropdown(
            input_plugin="quantumespresso.dos",
            description="dos.x:",
            setup_code_params={
                "computer": "localhost",
                "description": "dos.x in AiiDAlab container.",
                "label": "dos-6.7",
                "input_plugin": "quantumespresso.dos",
                "remote_abs_path": "/usr/bin/dos.x",
            },
        )
        self.projwfc = CodeDropdown(
            input_plugin="quantumespresso.projwfc",
            description="projwfc.x:",
            setup_code_params={
                "computer": "localhost",
                "description": "projwfc.x in AiiDAlab container.",
                "label": "projwfc-6.7",
                "input_plugin": "quantumespresso.projwfc",
                "remote_abs_path": "/usr/bin/projwfc.x",
            },
        )
        super().__init__(
            children=[
                self.codes_title,
                self.codes_help,
                self.pw,
                self.dos,
                self.projwfc,
            ],
            **kwargs,
        )


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    # This number provides a rough estimate for how many MPI tasks are needed
    # for a given structure.
    NUM_SITES_PER_MPI_TASK_DEFAULT = 6

    # Warn the user if they are trying to run calculations for a large
    # structure on localhost.
    RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD = 10

    input_structure = traitlets.Instance(StructureData, allow_none=True)
    process = traitlets.Instance(WorkChainNode, allow_none=True)
    disabled = traitlets.Bool()
    expert_mode = traitlets.Bool()
    _submission_blockers = traitlets.List(traitlets.Unicode)

    def __init__(self, **kwargs):
        self.message_area = ipw.Output()
        self._submission_blocker_messages = ipw.HTML()

        self.workchain_settings = WorkChainSettings()
        self.workchain_settings.relax_type.observe(self._update_state, "value")
        self.workchain_settings.bands_run.observe(self._update_state, "value")
        self.workchain_settings.pdos_run.observe(self._update_state, "value")

        self.kpoints_settings = KpointSettings()
        self.pseudo_family_selector = PseudoFamilySelector()
        self.codes_selector = CodeSettings()
        self.resources_config = ResourceSelectionWidget()

        self.set_input_parameters(DEFAULT_PARAMETERS)

        self.codes_selector.pw.observe(self._update_state, "selected_code")
        self.codes_selector.pw.observe(
            self._set_num_mpi_tasks_to_default, "selected_code"
        )
        self.codes_selector.dos.observe(self._update_state, "selected_code")
        self.codes_selector.projwfc.observe(self._update_state, "selected_code")
        ipw.dlink(
            (self.workchain_settings.workchain_protocol, "value"),
            (self.kpoints_settings, "kpoints_distance_default"),
            lambda protocol: PwBaseWorkChain.get_protocol_inputs(protocol)[
                "kpoints_distance"
            ],
        )

        self.tab = ipw.Tab(
            children=[
                self.workchain_settings,
            ],
            layout=ipw.Layout(min_height="250px"),
        )

        self.tab.set_title(0, "Workflow")

        self.submit_button = ipw.Button(
            description="Submit",
            tooltip="Submit the calculation with the selected parameters.",
            icon="play",
            button_style="success",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )

        self.submit_button.on_click(self._on_submit_button_clicked)

        self.expert_mode_control = ipw.ToggleButton(
            description="Expert mode",
            tooltip="Activate Expert mode for access to advanced settings.",
            value=True,
        )
        ipw.link((self, "expert_mode"), (self.expert_mode_control, "value"))

        # The SSSP installation status widget shows the installation status of
        # the SSSP pseudo potentials and triggers the installation in case that
        # they are not yet installed. The widget will remain in a "busy" state
        # in case that the installation was already triggered elsewhere, e.g.,
        # by the start up scripts.  The submission is blocked while the
        # potentials are not yet installed.
        self.sssp_installation_status = SSSPInstallWidget()
        self.sssp_installation_status.observe(self._update_state, ["busy", "installed"])
        self.sssp_installation_status.observe(self._toggle_install_widgets, "installed")

        # The QE setup widget checks whether there are codes that match specific
        # expected labels (e.g. "pw-6.7@localhost") and triggers both the
        # installation of QE into a dedicated conda environment and the setup of
        # the codes in case that they are not already configured.
        self.qe_setup_status = QESetupWidget()
        self.qe_setup_status.observe(self._update_state, "busy")
        self.qe_setup_status.observe(self._toggle_install_widgets, "installed")
        self.qe_setup_status.observe(self._auto_select_code, "installed")

        super().__init__(
            children=[
                self.message_area,
                self.tab,
                self.sssp_installation_status,
                self.qe_setup_status,
                self._submission_blocker_messages,
                ipw.HBox([self.submit_button, self.expert_mode_control]),
            ]
        )

    @traitlets.observe("expert_mode")
    def _observe_expert_mode(self, change):
        if change["new"]:
            self.tab.set_title(0, "Workflow")
            self.tab.set_title(1, "Advanced settings")
            self.tab.set_title(2, "Codes & Resources")
            self.tab.children = [
                self.workchain_settings,
                ipw.VBox(children=[self.pseudo_family_selector, self.kpoints_settings]),
                ipw.VBox(children=[self.codes_selector, self.resources_config]),
            ]
        else:
            self.tab.set_title(0, "Workflow")
            self.tab.children = [
                self.workchain_settings,
            ]

    @traitlets.observe("_submission_blockers")
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
        # No input structure specified.
        if self.input_structure is None:
            yield "No structure selected."

        # Do not submit while any of the background setup processes are running.
        if self.qe_setup_status.busy or self.sssp_installation_status.busy:
            yield "Background setup processes must finish."

        # No code selected (this is ignored while the setup process is running).
        if (
            self.codes_selector.pw.selected_code is None
            and not self.qe_setup_status.busy
        ):
            yield (
                'No pw code selected. Go to "Expert mode" and then "Codes & '
                'Resources" to select a code.'
            )

        # No code selected for pdos (this is ignored while the setup process is running).
        if (
            self.workchain_settings.pdos_run.value
            and (
                self.codes_selector.dos.selected_code is None
                or self.codes_selector.projwfc.selected_code is None
            )
            and not self.qe_setup_status.busy
        ):
            yield (
                'No code selected to run PDOS. Go to "Expert mode" and then "Codes & '
                'Resources" to select a code.'
            )

        # SSSP library not installed
        if not self.sssp_installation_status.installed:
            yield "The SSSP library is not installed."

        # Would not actually run anything
        if not (
            self.workchain_settings.relax_type.value != "none"
            or self.workchain_settings.bands_run.value
            or self.workchain_settings.pdos_run.value
        ):
            yield (
                "Select either a structure relaxation method or at least one of the "
                "the bands or the PDOS calculations or both."
            )

    def _get_state(self):
        # Process is already running.
        if self.process is not None:
            return self.State.SUCCESS

        # Input structure not specified.
        if self.input_structure is None:
            self._submission_blockers = ["No structure selected."]
            # This blocker is handled differently than the other blockers,
            # because it is displayed as INIT state.
            return self.State.INIT

        blockers = list(self._identify_submission_blockers())
        if any(blockers):
            self._submission_blockers = blockers
            return self.State.READY
        else:
            self._submission_blockers = []
            return self.state.CONFIGURED

    def _update_state(self, _=None):
        self.state = self._get_state()

    def _toggle_install_widgets(self, change):
        if change["new"]:
            self.children = [
                child for child in self.children if child is not change["owner"]
            ]

    def _auto_select_code(self, change):
        if change["new"] and not change["old"]:
            for selector, code in [
                ("pw", "pw_code"),
                ("dos", "dos_code"),
                ("projwfc", "projwfc_code"),
            ]:
                try:
                    code_widget = getattr(self.codes_selector, selector)
                    code_widget.refresh()
                    code_widget.selected_code = load_code(DEFAULT_PARAMETERS[code])
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

    def _get_default_num_mpi_tasks(self):
        """Determine a reasonable value for the number of MPI tasks for the selected structure."""
        if self.codes_selector.pw.selected_code:
            num_sites = len(self.input_structure.sites) if self.input_structure else 1
            num_mpi_tasks = max(
                1, ceil(num_sites / self.NUM_SITES_PER_MPI_TASK_DEFAULT)
            )
            return num_mpi_tasks

        return 1

    def _set_num_mpi_tasks_to_default(self, _=None):
        """Set the number of MPI tasks to a reasonable value for the selected structure."""
        self.resources_config.num_mpi_tasks.value = self._get_default_num_mpi_tasks()
        self._check_resources()

    def _check_resources(self):
        """Check whether the currently selected resources will be sufficient and warn if not."""
        if not self.codes_selector.pw.selected_code:
            return  # No code selected, nothing to do.

        num_mpi_tasks = self.resources_config.num_mpi_tasks.value
        on_localhost = (
            self.codes_selector.pw.selected_code.computer.get_hostname() == "localhost"
        )
        if self.codes_selector.pw.selected_code and on_localhost and num_mpi_tasks > 1:
            self._show_alert_message(
                "The selected code would be executed on the local host, but "
                "the number of MPI tasks is larger than one. Please review "
                "the configuration and consider to select a code that runs "
                'on a larger system if necessary (see the "Codes & '
                'Resources" tab).',
                alert_class="warning",
            )
            self.expert_mode = True
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
                'necessary (see the "Codes & Resources" tab).',
                alert_class="warning",
            )
            self.expert_mode = True

    def _get_cpus_per_node(self):
        """Determine the default number of CPUs per node based on the code configuration."""
        if self.codes_selector.pw.selected_code:
            selected_code = self.codes_selector.pw.selected_code
            return selected_code.computer.metadata["default_mpiprocs_per_machine"]
        return 1

    def _determine_resources(self):
        """Calculate the number of nodes and tasks per node."""
        cpus_per_node = self._get_cpus_per_node()
        num_mpi_tasks_selected = self.resources_config.num_mpi_tasks.value

        num_nodes = max(1, ceil(num_mpi_tasks_selected / cpus_per_node))
        num_mpi_tasks_per_node = ceil(num_mpi_tasks_selected / num_nodes)

        return {
            "num_machines": num_nodes,
            "num_mpiprocs_per_machine": num_mpi_tasks_per_node,
        }

    @traitlets.observe("state")
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.disabled = change["new"] not in (
                self.State.READY,
                self.State.CONFIGURED,
            )
            self.submit_button.disabled = change["new"] != self.State.CONFIGURED

    @traitlets.observe("input_structure")
    def _observe_input_structure(self, change):
        self.set_input_parameters(DEFAULT_PARAMETERS)
        self._update_state()
        self._set_num_mpi_tasks_to_default()

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                self.input_structure = process_node.inputs.structure
                builder_parameters = process_node.get_extra("builder_parameters", None)
                if builder_parameters is not None:
                    self.set_input_parameters(builder_parameters)
            self._update_state()

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.submit()

    def get_input_parameters(self):
        """Get the builder parameters based on the GUI inputs."""

        parameters = dict(
            # Work chain settings
            relax_type=self.workchain_settings.relax_type.value,
            electronic_type=self.workchain_settings.electronic_type.value,
            spin_type=self.workchain_settings.spin_type.value,
            run_bands=self.workchain_settings.bands_run.value,
            run_pdos=self.workchain_settings.pdos_run.value,
            protocol=self.workchain_settings.workchain_protocol.value,
            # Codes
            pw_code=self.codes_selector.pw.selected_code.uuid,
            dos_code=self.codes_selector.dos.selected_code.uuid,
            projwfc_code=self.codes_selector.projwfc.selected_code.uuid,
            # Advanced settings
            pseudo_family=self.pseudo_family_selector.value,
        )
        if self.kpoints_settings.override_protocol_kpoints.value:
            parameters[
                "kpoints_distance_override"
            ] = self.kpoints_settings.kpoints_distance.value

        return parameters

    def set_input_parameters(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        # Codes
        def _load_code(code):
            if code is not None:
                try:
                    return load_code(code)
                except NotExistent:
                    return None

        with self.hold_trait_notifications():
            # Work chain settings
            self.workchain_settings.relax_type.value = parameters["relax_type"]
            self.workchain_settings.spin_type.value = parameters["spin_type"]
            self.workchain_settings.electronic_type.value = parameters[
                "electronic_type"
            ]
            self.workchain_settings.bands_run.value = parameters["run_bands"]
            self.workchain_settings.pdos_run.value = parameters["run_pdos"]
            self.workchain_settings.workchain_protocol.value = parameters["protocol"]
            # Codes
            self.codes_selector.pw.selected_code = _load_code(parameters["pw_code"])
            if parameters["run_pdos"]:
                self.codes_selector.dos.selected_code = _load_code(
                    parameters["dos_code"]
                )
                self.codes_selector.projwfc.selected_code = _load_code(
                    parameters["projwfc_code"]
                )
            # Advanced settings
            self.pseudo_family_selector.value = parameters["pseudo_family"]
            if parameters.get("kpoints_distance_override", None) is not None:
                self.kpoints_settings.kpoints_distance.value = parameters[
                    "kpoints_distance_override"
                ]
                self.kpoints_settings.override_protocol_kpoints.value = True

    def submit(self, _=None):
        assert self.input_structure is not None
        parameters = self.get_input_parameters()

        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            pw_code=load_code(parameters["pw_code"]),
            dos_code=load_code(parameters["dos_code"]),
            projwfc_code=load_code(parameters["projwfc_code"]),
            protocol=parameters["protocol"],
            pseudo_family=parameters["pseudo_family"],
            relax_type=RelaxType(parameters["relax_type"]),
            spin_type=SpinType(parameters["spin_type"]),
            electronic_type=ElectronicType(parameters["electronic_type"]),
        )

        if "kpoints_distance_override" in parameters:
            builder.kpoints_distance_override = Float(
                parameters["kpoints_distance_override"]
            )

        if not parameters.pop("run_bands"):
            builder.pop("bands")

        if not parameters.pop("run_pdos"):
            builder.pop("pdos")

        resources = self._determine_resources()
        update_resources(builder, resources)

        self.process = submit(builder)

        # Set the builder parameters on the work chain
        self.process.set_extra("builder_parameters", self.get_input_parameters())

    def reset(self):
        with self.hold_trait_notifications():
            self.process = None
            self.input_structure = None
            self.set_input_parameters(DEFAULT_PARAMETERS)


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppWidgetStep):

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        self.process_tree = ProcessNodesTreeWidget()
        ipw.dlink((self, "process"), (self.process_tree, "process"))

        self.node_view = NodeViewWidget(layout={"width": "auto", "height": "auto"})
        ipw.dlink(
            (self.process_tree, "selected_nodes"),
            (self.node_view, "node"),
            transform=lambda nodes: nodes[0] if nodes else None,
        )
        self.process_status = ipw.VBox(children=[self.process_tree, self.node_view])

        # Setup process monitor
        self.process_monitor = ProcessMonitor(
            timeout=0.2,
            callbacks=[
                self.process_tree.update,
                self._update_state,
            ],
        )
        ipw.dlink((self, "process"), (self.process_monitor, "process"))

        super().__init__([self.process_status], **kwargs)

    def can_reset(self):
        "Do not allow reset while process is running."
        return self.state is not self.State.ACTIVE

    def reset(self):
        self.process = None

    def _update_state(self):
        if self.process is None:
            self.state = self.State.INIT
        else:
            process_state = self.process.process_state
            if process_state in (
                ProcessState.CREATED,
                ProcessState.RUNNING,
                ProcessState.WAITING,
            ):
                self.state = self.State.ACTIVE
            elif process_state in (ProcessState.EXCEPTED, ProcessState.KILLED):
                self.state = self.State.FAIL
            elif process_state is ProcessState.FINISHED:
                self.state = self.State.SUCCESS

    @traitlets.observe("process")
    def _observe_process(self, change):
        self._update_state()
