"""Widgets for the submission of bands work chains.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import ipywidgets as ipw
import traitlets
from aiida.engine import ProcessState
from aiida.engine import submit
from aiidalab_widgets_base import CodeDropdown
from aiida.orm import ProcessNode
from aiida.plugins import DataFactory

from process import ProcessMonitor
from process import ProcessNodesTreeWidget
from pseudos import PseudoFamilySelector
from widgets import NodeViewWidget
from widgets import ResourceSelectionWidget
from wizard import WizardApp, WizardAppStep

from IPython.display import clear_output, display

from apps.quantumespresso.qe_workflow import QeAppWorkChain
from aiida_quantumespresso.common.types import SpinType, ElectronicType, RelaxType

StructureData = DataFactory("structure")

WARNING_ICON = "\u26A0"


def update_resources(builder, resources):
    for k, v in builder.items():
        if isinstance(v, dict):
            if k == "resources":
                builder["resources"].update(resources)
            else:
                update_resources(v, resources)


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppStep):
    """Step for submission of a bands workchain."""

    input_structure = traitlets.Instance(StructureData, allow_none=True)

    process = traitlets.Instance(ProcessNode, allow_none=True)
    disabled = traitlets.Bool()

    def __init__(self, description=None, pseudo=None, kpoints_distance=None, **kwargs):
        self.code_group_pw = CodeDropdown(
            input_plugin="quantumespresso.pw",
            text="Select PW code",
            setup_code_params={
                "computer": "localhost",
                "description": "pw.x in AiiDAlab container.",
                "label": "pw",
                "input_plugin": "quantumespresso.pw",
                "remote_abs_path": "/usr/bin/pw.x",
            },
        )

        self.code_group_pw.observe(lambda _: self._update_state(), ["selected_code"])

        # Setup the compute resources tab
        self.resources = ResourceSelectionWidget()

        # Clicking on the 'submit' button will trigger the execution of the
        # submit() method.
        self.submit_button = ipw.Button(
            description="Submit",
            tooltip="Submit the calculation with the selected parameters.",
            icon="play",
            button_style="success",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )
        self.submit_button.on_click(self._on_submit_button_clicked)

        # Setup pseudo potential family selection
        self.pseudo_family_selector = PseudoFamilySelector()

        # Show warning in cofig title when pseudos are not installed:
        pseudo_output = ipw.Output()

        def _observe_pseudo_modify(value):
            with pseudo_output:
                clear_output()
                if value["new"]:
                    display(self.pseudo_family_selector)

        self.modify_pseudo = ipw.Checkbox(
            description="Non-default pseudo",
            indent=False,
        )
        self.modify_pseudo.observe(_observe_pseudo_modify, "value")

        if pseudo:
            self.modify_pseudo.value = True
            self.pseudo_family_selector.value = pseudo

        def _observe_sssp_installed(change):
            self._observe_state(change=dict(new=self.state))  # trigger refresh

        self.pseudo_family_selector.observe(_observe_sssp_installed, "installed")
        _observe_sssp_installed(
            change=dict(new=self.pseudo_family_selector.installed)
        )  # init

        ipw.dlink((self, "disabled"), (self.code_group_pw.dropdown, "disabled"))
        ipw.dlink((self, "disabled"), (self.resources.number_of_nodes, "disabled"))
        ipw.dlink((self, "disabled"), (self.resources.cpus_per_node, "disabled"))
        ipw.dlink((self, "disabled"), (self.pseudo_family_selector, "disabled"))

        # Initialize widget disabled status based on step state.
        self.disabled = self.state != WizardApp.State.READY

        # Modify k-points distance.
        self.kpoints_distance = ipw.FloatText(
            value=0.5, description="K-points distance:", disabled=False
        )
        kpoints_output = ipw.Output()

        def _observe_kpoints_modify(value):
            with kpoints_output:
                clear_output()
                if value["new"]:
                    display(self.kpoints_distance)

        self.modify_kpoints_distance = ipw.Checkbox(
            description="Non-default k-points distance",
            indent=False,
        )
        self.modify_kpoints_distance.observe(_observe_kpoints_modify, "value")

        if kpoints_distance:
            self.modify_kpoints_distance.value = True

        # Spin type.
        self.spin_type = ipw.Dropdown(
            options=[t.name for t in SpinType],
            value="NONE",
            description="Spin Type:",
        )

        # Electronic type.
        self.electronic_type = ipw.Dropdown(
            options=[t.name for t in ElectronicType],
            value="METAL",
            description="Electronic Type:",
        )

        # Simulation protocol.
        self.simulation_protocol = ipw.Dropdown(
            options=["fast", "moderate", "precise"],
            value="moderate",
            description="Protocol",
        )

        # Geometry optimization.
        self.run_geo_opt = ipw.Dropdown(
            options=[t.name for t in RelaxType],
            value="POSITIONS",
            description="Geometry Optimization:",
        )

        # Band Structure.
        self.run_bands = ipw.Checkbox(
            value=True, description="Run Band Structure", disabled=False, indent=False
        )

        # PDOS.
        self.run_pdos = ipw.Checkbox(
            value=False,
            description="Run PDOS",
            indent=False,
        )

        self.code_group_dos = CodeDropdown(
            input_plugin="quantumespresso.dos",
            text="Select DOS code",
            setup_code_params={
                "computer": "localhost",
                "description": "dos.x in AiiDAlab container.",
                "label": "dos",
                "input_plugin": "quantumespresso.dos",
                "remote_abs_path": "/usr/bin/dos.x",
            },
        )

        self.code_group_projwfc = CodeDropdown(
            input_plugin="quantumespresso.projwfc",
            text="Select PROJWFC code",
            setup_code_params={
                "computer": "localhost",
                "description": "projwfc.x in AiiDAlab container.",
                "label": "projwfc",
                "input_plugin": "quantumespresso.projwfc",
                "remote_abs_path": "/usr/bin/projwfc.x",
            },
        )

        pdos_code_output = ipw.Output()

        def _code_setup_visibility(value):
            with pdos_code_output:
                clear_output()
                if value["new"]:
                    display(self.code_group_dos, self.code_group_projwfc)

        self.run_pdos.observe(_code_setup_visibility, names="value")

        super().__init__(
            children=[
                ipw.Label(
                    'Specify the parameters and options for the calculation and then click on "Submit".'
                ),
                self.code_group_pw,
                self.modify_pseudo,
                pseudo_output,
                self.modify_kpoints_distance,
                kpoints_output,
                self.resources,
                self.simulation_protocol,
                self.electronic_type,
                self.spin_type,
                self.run_geo_opt,
                self.run_bands,
                self.run_pdos,
                pdos_code_output,
                self.submit_button,
            ],
            **kwargs,
        )

    def _update_state(self):
        if self.process is None:
            if self.input_structure is None:
                self.state = WizardApp.State.INIT
            elif (
                self.code_group_pw.selected_code is None
                or not self.pseudo_family_selector.installed
            ):
                self.state = WizardApp.State.READY
            else:
                self.state = WizardApp.State.CONFIGURED
        else:
            self.state = WizardApp.State.SUCCESS

    @traitlets.observe("state")
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.disabled = change["new"] not in (
                WizardApp.State.READY,
                WizardApp.State.CONFIGURED,
            )
            self.submit_button.disabled = change["new"] != WizardApp.State.CONFIGURED

    @traitlets.observe("input_structure")
    def _observe_input_structure(self, change):
        self._update_state()

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                self.input_structure = process_node.inputs.structure

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.state = WizardApp.State.ACTIVE
        self.submit()

    def submit(self, _=None):
        assert self.input_structure is not None

        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            pw_code=self.code_group_pw.selected_code,
            dos_code=self.code_group_dos.selected_code,
            projwfc_code=self.code_group_projwfc.selected_code,
            protocol=self.simulation_protocol.value,
            relax_type=RelaxType[self.run_geo_opt.value],
            spin_type=SpinType[self.spin_type.value],
            electronic_type=ElectronicType[self.electronic_type.value],
            pseudo_family=self.pseudo_family_selector.value
            if self.modify_pseudo.value
            else None,
            kpoints_distance_override=self.kpoints_distance
            if self.modify_kpoints_distance
            else None,
        )

        if not self.run_bands.value:
            builder.pop("bands")
        if not self.run_pdos.value:
            builder.pop("pdos")

        resources = {
            "num_machines": self.resources.number_of_nodes.value,
            "num_mpiprocs_per_machine": self.resources.cpus_per_node.value,
        }
        update_resources(builder, resources)

        self.process = submit(builder)

    def reset(self):
        self.process = None


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppStep):

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        self.process_tree = ProcessNodesTreeWidget(
            refresh_period=-1,  # managed within this class,
        )
        ipw.dlink((self, "process"), (self.process_tree, "process"))

        self.node_view = NodeViewWidget(
            layout={"width": "auto", "height": "auto", "border": "1px solid black"}
        )
        ipw.dlink(
            (self.process_tree, "selected_nodes"),
            (self.node_view, "node"),
            transform=lambda nodes: nodes[0] if nodes else None,
        )
        self.process_status = ipw.VBox(children=[self.process_tree, self.node_view])

        # Setup process monitor
        self.process_monitor = ProcessMonitor(
            timeout=0.2,  # run every half second
            callbacks=[
                (self.process_tree.update, 1),
            ],
        )
        ipw.dlink((self, "process"), (self.process_monitor, "process"))

        super().__init__([self.process_status], **kwargs)

    def reset(self):
        self.process = None

    def _update_state(self):
        if self.process is None:
            self.state = WizardApp.State.INIT
        else:
            process_state = self.process.process_state
            if process_state in (
                ProcessState.CREATED,
                ProcessState.RUNNING,
                ProcessState.WAITING,
            ):
                self.state = WizardApp.State.ACTIVE
            elif process_state in (ProcessState.EXCEPTED, ProcessState.KILLED):
                self.state = WizardApp.State.FAIL
            elif process_state is ProcessState.FINISHED:
                self.state = WizardApp.State.SUCCESS

    @traitlets.observe("process")
    def _observe_process(self, change):
        self._update_state()
