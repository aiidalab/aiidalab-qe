"""Widgets for the submission of bands work chains.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
from math import ceil
from pprint import pformat

import ipywidgets as ipw
import traitlets
from aiida.common import NotExistent
from aiida.engine import ProcessState, submit
from aiida.orm import ProcessNode, WorkChainNode, load_code
from aiida.plugins import DataFactory
from aiida_quantumespresso.common.types import ElectronicType, RelaxType, SpinType
from aiidalab_widgets_base import (
    CodeDropdown,
    ProcessMonitor,
    ProcessNodesTreeWidget,
    WizardAppWidgetStep,
)

from aiidalab_qe.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.pseudos import PseudoFamilySelector
from aiidalab_qe.widgets import NodeViewWidget, ResourceSelectionWidget
from aiidalab_qe_workchain import QeAppWorkChain

StructureData = DataFactory("structure")


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
    button_style_on = "info"
    button_style_off = "danger"

    spin_type = traitlets.Instance(SpinType, allow_none=True)
    electronic_type = traitlets.Instance(ElectronicType, allow_none=True)

    def __init__(self, **kwargs):

        # Work chain protocol.
        self.geo_opt_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "NONE"),
                ("Atomic positions", "POSITIONS"),
                ("Full geometry", "POSITIONS_CELL"),
            ],
            value="POSITIONS_CELL",
        )
        self._spin_type = ipw.Dropdown(
            options=[("Non-magnetic", "NONE"), ("Ferromagnetic", "COLLINEAR")],
            value=DEFAULT_PARAMETERS["spin_type"].upper(),
            description="Magnetism:",
            style={"description_width": "initial"},
        )
        self._spin_type.observe(self.set_spin_type_trait, "value")
        self.set_spin_type_trait()

        self._electronic_type = ipw.Dropdown(
            options=[("Metal", "METAL"), ("Insulator", "INSULATOR")],
            value=DEFAULT_PARAMETERS["electronic_type"].upper(),
            description="Electronic Type:",
            style={"description_width": "initial"},
        )
        self._electronic_type.observe(self.set_electronic_type_trait, "value")
        self.set_electronic_type_trait()

        self.bands_run = ipw.Checkbox(
            description="",
            tooltip="Calculate the electronic band structure.",
            indent=False,
            value=True,
            layout=ipw.Layout(max_width="10%"),
        )

        # Work chain protocol.
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        super().__init__(
            children=[
                self.structure_title,
                self.structure_help,
                self.geo_opt_type,
                self.materials_help,
                self._spin_type,
                self._electronic_type,
                self.properties_title,
                ipw.HTML("Select which properties to calculate:"),
                ipw.HBox(children=[ipw.HTML("<b>Band structure</b>"), self.bands_run]),
                self.properties_help,
                self.protocol_title,
                ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
                self.workchain_protocol,
                self.protocol_help,
            ],
            **kwargs,
        )

    def set_spin_type_trait(self, _=None):
        self.spin_type = SpinType[self._spin_type.value.upper()]

    def set_electronic_type_trait(self, _=None):
        self.electronic_type = ElectronicType[self._electronic_type.value]


class KpointSettings(ipw.VBox):

    # The default of `kpoints_distance` must be linked to the `protocol`
    kpoints_distance_default = traitlets.Float(default_value=0.15)
    kpoints_distance = traitlets.Float(allow_none=True)
    degauss = traitlets.Float(allow_none=True)

    kpoints_distance_description = ipw.HTML(
        """<p>
        The k-points mesh density is set by the chosen <b>protocol</b>.
        Untick the box to override the default.
    </p>"""
    )

    def __init__(self, **kwargs):

        self._set_kpoints_distance_automatically = ipw.Checkbox(
            description="Use default k-points distance.",
            indent=False,
            value=True,
        )
        self._kpoints_distance = ipw.FloatText(
            value=self.kpoints_distance_default,
            step=0.05,
            description="K-points distance:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._set_kpoints_distance_automatically, "value"),
            (self._kpoints_distance, "disabled"),
        )

        self._kpoints_distance.observe(self.set_kpoints_distance_trait, "value")
        self._set_kpoints_distance_automatically.observe(
            self.set_kpoints_distance_trait, "value"
        )
        self.set_kpoints_distance_trait()

        super().__init__(
            children=[
                self.kpoints_distance_description,
                ipw.HBox(
                    [self._set_kpoints_distance_automatically, self._kpoints_distance]
                ),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    def set_kpoints_distance_trait(self, _=None):
        self.kpoints_distance = (
            self.kpoints_distance_default
            if self._set_kpoints_distance_automatically.value
            else self._kpoints_distance.value
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
                "label": "pw",
                "input_plugin": "quantumespresso.pw",
                "remote_abs_path": "/usr/bin/pw.x",
            },
        )
        super().__init__(
            children=[
                self.codes_title,
                self.codes_help,
                self.pw,
            ],
            **kwargs,
        )


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    # The app will issue a warning to the user if the ratio between the total
    # number of sites and the total number of CPUs is larger than this value:
    NUM_SITES_PER_MPI_TASK_DEFAULT = 6

    input_structure = traitlets.Instance(StructureData, allow_none=True)
    process = traitlets.Instance(WorkChainNode, allow_none=True)
    disabled = traitlets.Bool()
    builder_parameters = traitlets.Dict()
    expert_mode = traitlets.Bool()

    def __init__(self, **kwargs):
        self.message_area = ipw.Output()
        self.workchain_settings = WorkChainSettings()
        self.kpoints_settings = KpointSettings()
        self.pseudo_family_selector = PseudoFamilySelector()
        self.codes_selector = CodeSettings()
        self.resources_config = ResourceSelectionWidget()

        self.set_trait("builder_parameters", self._default_builder_parameters())
        self._setup_builder_parameters_update()

        self.codes_selector.pw.observe(self._update_state, "selected_code")
        self.codes_selector.pw.observe(
            self._set_num_mpi_tasks_to_default, "selected_code"
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

        self._update_builder_parameters()

        self.builder_parameters_view = ipw.HTML(layout=ipw.Layout(width="auto"))
        ipw.dlink(
            (self, "builder_parameters"),
            (self.builder_parameters_view, "value"),
            transform=lambda p: '<pre style="line-height: 100%">'
            + pformat(p, indent=2, width=200)
            + "</pre>",
        )

        super().__init__(
            children=[
                self.message_area,
                self.tab,
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

    def _get_state(self):

        # Process is already running.
        if self.process is not None:
            return self.State.SUCCESS

        # Input structure not specified.
        if self.input_structure is None:
            return self.State.INIT

        # PW code not selected.
        if self.codes_selector.pw.selected_code is None:
            return self.State.READY

        return self.State.CONFIGURED

    def _update_state(self, _=None):
        self.state = self._get_state()

    _ALERT_MESSAGE = """
        <div class="alert alert-{alert_class} alert-dismissible">
        <a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>
        <span class="closebtn" onclick="this.parentElement.style.display='none';">&times;</span>
        <strong>{message}</strong>
        </div>"""

    def _show_alert_message(self, message, alert_class="info"):
        with self.message_area:
            display(  # noqa
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
        num_mpi_tasks = self.resources_config.num_mpi_tasks.value
        if self.codes_selector.pw.selected_code and num_mpi_tasks > 1:
            hostname = self.codes_selector.pw.selected_code.computer.get_hostname()
            if hostname == "localhost":
                self._show_alert_message(
                    "The selected code would be executed on the localhost, but "
                    "the number of MPI tasks is larger than one. Please review "
                    "the configuration and consider to select a code that runs "
                    'on a larger system if necessary (see the "Codes & '
                    'Resources" tab).',
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
        self.set_trait("builder_parameters", self._default_builder_parameters())
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
                    self.set_trait("builder_parameters", builder_parameters)
            self._update_state()

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.submit()

    def _setup_builder_parameters_update(self):
        """Set up all ``observe`` calls to monitor changes in user inputs."""
        update = self._update_builder_parameters  # alias for code conciseness
        # Properties
        self.workchain_settings.geo_opt_type.observe(update, ["value"])
        self.workchain_settings.bands_run.observe(update, ["value"])
        # Codes
        self.codes_selector.pw.observe(update, ["selected_code"])
        # Material settings
        self.workchain_settings.observe(update, ["electronic_type"])
        self.workchain_settings.observe(update, ["spin_type"])
        # Calculation settings
        self.workchain_settings.workchain_protocol.observe(update, ["value"])
        self.kpoints_settings.observe(update, ["kpoints_distance"])
        self.pseudo_family_selector.observe(update, ["value"])

    @staticmethod
    def _serialize_builder_parameters(parameters):
        parameters = parameters.copy()  # create copy to not modify original dict

        # Codes
        def _get_uuid(code):
            return None if code is None else str(code.uuid)

        parameters["pw_code"] = _get_uuid(parameters["pw_code"])

        # Protocol
        parameters["electronic_type"] = parameters["electronic_type"].value
        parameters["relax_type"] = parameters["relax_type"].value
        parameters["spin_type"] = parameters["spin_type"].value
        return parameters

    @staticmethod
    def _deserialize_builder_parameters(parameters):
        parameters = parameters.copy()  # create copy to not modify original dict

        # Codes
        def _load_code(code):
            if code is not None:
                try:
                    return load_code(code)
                except NotExistent as error:
                    print("error", error)
                    return None

        parameters["pw_code"] = _load_code(parameters["pw_code"])

        # Protocol
        parameters["electronic_type"] = ElectronicType(parameters["electronic_type"])
        parameters["relax_type"] = RelaxType(parameters["relax_type"])
        parameters["spin_type"] = SpinType(parameters["spin_type"])
        return parameters

    def _update_builder_parameters(self, _=None):
        self.set_trait(
            "builder_parameters",
            self._serialize_builder_parameters(
                dict(
                    # Properties
                    relax_type=RelaxType[self.workchain_settings.geo_opt_type.value],
                    run_bands=self.workchain_settings.bands_run.value,
                    # Codes
                    pw_code=self.codes_selector.pw.selected_code,
                    # Material settings
                    electronic_type=self.workchain_settings.electronic_type,
                    spin_type=self.workchain_settings.spin_type,
                    # Calculation settings
                    kpoints_distance_override=self.kpoints_settings.kpoints_distance,
                    protocol=self.workchain_settings.workchain_protocol.value,
                    pseudo_family=self.pseudo_family_selector.value,
                )
            ),
        )

    @traitlets.observe("builder_parameters")
    def _observe_builder_parameters(self, change):
        bp = self._deserialize_builder_parameters(change["new"])

        with self.hold_trait_notifications():
            # Properties
            relax_type = bp.get("relax_type", RelaxType["NONE"])
            self.workchain_settings.geo_opt_type.value = relax_type.value.upper()
            self.workchain_settings.bands_run.value = bp["run_bands"]
            # Codes
            self.codes_selector.pw.selected_code = bp.get("pw_code")
            # Material settings
            self.workchain_settings.spin_type = bp.get("spin_type", SpinType["NONE"])
            self.workchain_settings.electronic_type = bp.get(
                "electronic_type", ElectronicType["METAL"]
            )
            # Calculation settings
            self.kpoints_settings.kpoints_distance = bp.get("kpoints_distance")
            self.workchain_settings.workchain_protocol.value = bp["protocol"]
            self.pseudo_family_selector.value = bp["pseudo_family"]

    def submit(self, _=None):

        assert self.input_structure is not None

        builder_parameters = self.builder_parameters.copy()

        run_bands = builder_parameters.pop("run_bands")

        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            **self._deserialize_builder_parameters(builder_parameters),
        )

        if not run_bands:
            builder.pop("bands")

        resources = self._determine_resources()
        update_resources(builder, resources)

        self.process = submit(builder)
        self.process.set_extra("builder_parameters", self.builder_parameters.copy())

    def reset(self):
        with self.hold_trait_notifications():
            self.process = None
            self.input_structure = None
            self.builder_parameters = self._default_builder_parameters()

    @traitlets.default("builder_parameters")
    def _default_builder_parameters(self):
        return DEFAULT_PARAMETERS


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
