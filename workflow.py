"""Widgets for the submission of bands work chains.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import ipywidgets as ipw
import traitlets
from aiida.engine import ProcessState
from aiida.engine import submit
from aiida.orm import ProcessNode
from aiida.plugins import DataFactory
from aiidalab_widgets_base import CodeDropdown
from aiidalab_widgets_base import ProcessMonitor
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base import WizardAppWidgetStep

from pseudos import PseudoFamilySelector
from widgets import NodeViewWidget
from widgets import ResourceSelectionWidget

from IPython.display import clear_output, display

from apps.quantumespresso.qe_workflow import QeAppWorkChain
from aiida_quantumespresso.common.types import SpinType, ElectronicType, RelaxType


StructureData = DataFactory("structure")


def update_resources(builder, resources):
    for k, v in builder.items():
        if isinstance(v, dict):
            if k == "resources":
                builder["resources"].update(resources)
            else:
                update_resources(v, resources)


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    input_structure = traitlets.Instance(StructureData, allow_none=True)

    process = traitlets.Instance(ProcessNode, allow_none=True)
    disabled = traitlets.Bool()

    def __init__(
        self,
        description=None,
        pseudo_family=None,
        kpoints_distance=None,
        electronic_type=None,
        spin_type=None,
        degauss=None,
        **kwargs
    ):

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

        # What to compute.

        # Geometry optimization.
        self.run_geo_opt = ipw.Checkbox(
            value=True,
            description="Run Geometry Optimization",
            disabled=False,
            indent=False,
        )

        self.geo_opt_type = ipw.Dropdown(
            options=["POSITIONS", "POSITIONS_CELL"],
            value="POSITIONS",
            description="Geometry Optimization:",
            style={"description_width": "initial"},
        )

        def _geo_type_visibility(value):
            self.geo_opt_type.layout.visibility = (
                "visible" if value["new"] else "hidden"
            )

        self.run_geo_opt.observe(_geo_type_visibility, "value")

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

        # Simulation protocol.
        self.simulation_protocol = ipw.Dropdown(
            options=["fast", "moderate", "precise"],
            value="moderate",
            description="Protocol",
        )

        # DFT functional.
        self.dft_functional = ipw.Dropdown(
            options=["PBE", "PBEsol"],
            value="PBE",
            description="DFT functional",
            style={"description_width": "initial"},
        )

        # Spin type.
        self.spin_type = ipw.Dropdown(
            options=["NONE", "COLLINEAR", "NON_COLLINEAR"],
            value="NONE",
            description="Spin Type:",
        )

        # Electronic type.
        self.electronic_type = ipw.Dropdown(
            options=["METAL", "INSULATOR"],
            value="METAL",
            description="Electronic Type:",
            style={"description_width": "initial"},
        )

        # Pseudo potential family selection.
        self.modify_pseudo = ipw.Checkbox(
            description="Choose pseudo",
            indent=False,
        )

        pseudo_output = ipw.Output()

        def _observe_pseudo_modify(value):
            with pseudo_output:
                clear_output()
                if value["new"]:
                    display(self.pseudo_family_selector)

        self.modify_pseudo.observe(_observe_pseudo_modify, "value")

        self.pseudo_family_selector = PseudoFamilySelector()
        ipw.dlink(
            (self.dft_functional, "value"), (self.pseudo_family_selector, "functional")
        )

        def _observe_sssp_installed(change):
            self._observe_state(change=dict(new=self.state))  # trigger refresh

        # Show warning in cofig title when pseudos are not installed:
        self.pseudo_family_selector.observe(_observe_sssp_installed, "installed")
        _observe_sssp_installed(change=dict(new=self.pseudo_family_selector.installed))

        # Modify k-points distance.
        self.kpoints_distance = ipw.FloatText(
            value=0.5,
            step=0.1,
            description="K-points distance:",
            disabled=False,
            style={"description_width": "initial"},
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

        # Modify degauss.
        self.degauss = ipw.FloatText(
            value=0.01,
            step=0.01,
            description="Gaussian spreading",
            disabled=False,
            style={"description_width": "initial"},
        )
        degauss_output = ipw.Output()

        def _observe_degauss_modify(value):
            with degauss_output:
                clear_output()
                if value["new"]:
                    display(self.degauss)

        self.modify_degauss = ipw.Checkbox(
            description="Non-default Gaussian spreading",
            indent=False,
        )
        self.modify_degauss.observe(_observe_degauss_modify, "value")

        # Manage compute resources.
        self.resources = ResourceSelectionWidget()

        # Manage codes.
        self.code_group_pw = CodeDropdown(
            input_plugin="quantumespresso.pw",
            description="PW code:",
            setup_code_params={
                "computer": "localhost",
                "description": "pw.x in AiiDAlab container.",
                "label": "pw",
                "input_plugin": "quantumespresso.pw",
                "remote_abs_path": "/usr/bin/pw.x",
            },
        )

        self.code_group_dos = CodeDropdown(
            input_plugin="quantumespresso.dos",
            description="DOS code",
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
            description="PROJWFC code",
            setup_code_params={
                "computer": "localhost",
                "description": "projwfc.x in AiiDAlab container.",
                "label": "projwfc",
                "input_plugin": "quantumespresso.projwfc",
                "remote_abs_path": "/usr/bin/projwfc.x",
            },
        )

        self.code_group_pw.observe(self._update_state, "selected_code")
        self.code_group_dos.observe(self._update_state, "selected_code")
        self.code_group_projwfc.observe(self._update_state, "selected_code")
        self.run_pdos.observe(self._update_state, "value")

        code_output = ipw.Output()

        def _code_setup_visibility(value):
            with code_output:
                clear_output()
                display(self.code_group_pw)
                if value["new"]:
                    display(self.code_group_dos, self.code_group_projwfc)

        self.run_pdos.observe(_code_setup_visibility, names="value")
        _code_setup_visibility({"new": None})

        ipw.dlink((self, "disabled"), (self.code_group_pw.dropdown, "disabled"))
        ipw.dlink((self, "disabled"), (self.resources.number_of_nodes, "disabled"))
        ipw.dlink((self, "disabled"), (self.resources.cpus_per_node, "disabled"))
        ipw.dlink((self, "disabled"), (self.pseudo_family_selector, "disabled"))

        # Initialize widget disabled status based on step state.
        self.disabled = self.state != self.State.READY

        # Use input parameters.
        if pseudo_family:
            _, _, functional, protocol = pseudo_family.split("/")
            self.modify_pseudo.value = True
            self.dft_functional.value = functional
            self.pseudo_family_selector.protocol = protocol

        if kpoints_distance:
            self.modify_kpoints_distance.value = True
            self.kpoints_distance.value = kpoints_distance

        if electronic_type:
            self.electronic_type.value = electronic_type

        if spin_type:
            self.spin_type.value = spin_type

        if degauss:
            self.modify_degauss.value = True
            self.degauss.value = degauss

        super().__init__(
            children=[
                ipw.Label("Specify which calculations to run."),
                ipw.HBox([self.run_geo_opt, self.geo_opt_type]),
                self.run_bands,
                self.run_pdos,
                self.simulation_protocol,
                self.dft_functional,
                self.electronic_type,
                self.spin_type,
                self.modify_pseudo,
                pseudo_output,
                self.modify_kpoints_distance,
                kpoints_output,
                self.modify_degauss,
                degauss_output,
                self.resources,
                code_output,
                self.submit_button,
            ],
            **kwargs,
        )

    def _get_state(self):

        # Input structure not specified.
        if self.input_structure is None:
            return self.State.INIT

        # Process is already running.
        if self.process is not None:
            return self.State.SUCCESS

        # Pseudo family is not installed.
        if not self.pseudo_family_selector.installed:
            return self.State.READY

        # PW code not selected.
        if self.code_group_pw.selected_code is None:
            return self.State.READY

        # PDOS run requested, but codes are not specified.
        if self.run_pdos.value:
            if self.code_group_dos.selected_code is None:
                return self.State.READY
            if self.code_group_projwfc.selected_code is None:
                return self.State.READY

        return self.State.CONFIGURED

    def _update_state(self, _=None):
        self.state = self._get_state()

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
        self._update_state()

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                self.input_structure = process_node.inputs.structure

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.state = self.State.ACTIVE
        self.submit()

    def submit(self, _=None):
        assert self.input_structure is not None

        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            pw_code=self.code_group_pw.selected_code,
            dos_code=self.code_group_dos.selected_code,
            projwfc_code=self.code_group_projwfc.selected_code,
            protocol=self.simulation_protocol.value,
            relax_type=RelaxType[self.geo_opt_type.value]
            if self.run_geo_opt.value
            else RelaxType["NONE"],
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
        with self.hold_trait_notifications():
            self.pseudo_family_selector.reset()
            self.resources.reset()


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppWidgetStep):

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        self.process_tree = ProcessNodesTreeWidget()
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
