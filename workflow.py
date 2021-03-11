"""Widgets for the submission of bands work chains.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import yaml
import ipywidgets as ipw
import traitlets
from aiida.engine import ProcessState
from aiida.engine import submit
from aiidalab_widgets_base import CodeDropdown
from aiida.orm import ProcessNode
from aiida.orm import Dict, Float, Str
from aiida.plugins import DataFactory

from process import ProcessMonitor
from process import ProcessNodesTreeWidget
from pseudos import PseudoFamilySelector
from widgets import NodeViewWidget
from widgets import ResourceSelectionWidget
from wizard import WizardApp, WizardAppStep

from apps.quantumespresso.qe_workflow import QeAppWorkChain

StructureData = DataFactory("structure")

WARNING_ICON = "\u26A0"


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppStep):
    """Step for submission of a bands workchain."""

    input_structure = traitlets.Instance(StructureData, allow_none=True)

    process = traitlets.Instance(ProcessNode, allow_none=True)
    disabled = traitlets.Bool()

    def __init__(self, description=None, **kwargs):
        self.code_group_pw = CodeDropdown(
            input_plugin="quantumespresso.pw",
            text="Select code",
            setup_code_params={
                "computer": "localhost",
                "description": "pw.x in AiiDAlab container.",
                "label": "pw",
                "input_plugin": "quantumespresso.pw",
                "remote_abs_path": "/usr/bin/pw.x",
            },
        )
        self.code_group_pw.observe(lambda _: self._update_state(), ["selected_code"])

        # Setup pseudo potential family selection
        self.pseudo_family_selector = PseudoFamilySelector()

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

        # Show warning in cofig title when pseudos are not installed:
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

        # Geometry optimization.
        self.run_geo_opt = ipw.Checkbox(
            value=True,
            description="Run Geometry Optimization",
            disabled=False,
            indent=False,
        )
        # Band Structure.
        self.run_band_structure = ipw.Checkbox(
            value=True, description="Run Band Structure", disabled=False, indent=False
        )
        # PDOS.
        self.run_pdos = ipw.Checkbox(
            value=True, description="Run PDOS", disabled=False, indent=False
        )

        super().__init__(
            children=[
                ipw.Label(
                    'Specify the parameters and options for the calculation and then click on "Submit".'
                ),
                self.code_group_pw,
                self.pseudo_family_selector,
                self.resources,
                self.run_geo_opt,
                self.run_band_structure,
                self.run_pdos,
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

    @traitlets.validate("process")
    def _validate_process(self, proposal):
        # process_node = proposal["value"]
        # builder = process_node.get_builder_restart()
        """
        try:
            # Check that parameters are consistent with what we would expect for the app.
            # This ensures that both processes we create with the app and those that are
            # loaded are consistent.
            assert process_node.process_label == "PwBandsWorkChain"
            assert builder.scf.pseudo_family.value == builder.bands.pseudo_family.value
            assert (
                dict(load_default_parameters())
                == dict(builder.scf.pw["parameters"])
                == dict(builder.bands.pw["parameters"])
            )
            assert builder.scf.kpoints_distance == Float(0.8)
        except AssertionError as error:
            raise traitlets.TraitError(
                f"Unable to set process, parameters are inconsistent: {error}"
            )
        """
        return proposal["value"]

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                # Restore the code parameters
                # builder = process_node.get_builder_restart()
                """
                self.pseudo_family_selector.value = builder.relax.base[
                    "pseudo_family"
                ].value
                try:
                    self.code_group_pw.selected_code = builder.relax.base["pw"][
                        "code"
                    ].full_label
                except traitlets.TraitError:
                    # It's not considered a critical issue if the exact code cannot
                    # be retrieved, however we set it to "None" in the interface to
                    # indicate that to the user.
                    self.code_group_pw.selected_code = None
                """
                self.input_structure = process_node.inputs.structure

    @property
    def options(self):
        return {
            "max_wallclock_seconds": 3600 * 2,
            "resources": {
                "num_machines": self.resources.number_of_nodes.value,
                "num_mpiprocs_per_machine": self.resources.cpus_per_node.value,
            },
        }

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.state = WizardApp.State.ACTIVE
        self.submit()

    def submit(self, _=None):
        assert self.input_structure is not None

        builder = QeAppWorkChain.get_builder()

        # Input structure.
        builder.structure = self.input_structure

        # Geometry optimization section.
        if self.run_geo_opt.value:
            builder.relax.base.pw.code = self.code_group_pw.selected_code
            with open("parameters/relax.yaml") as file:
                builder.relax.base.pw.parameters = Dict(
                    dict=yaml.load(file, Loader=yaml.FullLoader)
                )
            builder.relax.base.pw.metadata.options = self.options
            builder.relax.base.kpoints_distance = Float(0.8)
            builder.relax.base.pseudo_family = Str(self.pseudo_family_selector.value)

        # Bands section.
        if self.run_band_structure.value:
            with open("parameters/bands.yaml") as file:
                builder.bands.scf.pw.parameters = Dict(
                    dict=yaml.load(file, Loader=yaml.FullLoader)
                )
            builder.bands.scf.pw.metadata.options = self.options
            builder.bands.scf.kpoints_distance = Float(0.8)
            builder.bands.scf.pseudo_family = Str(self.pseudo_family_selector.value)
            builder.bands.scf.pw.code = self.code_group_pw.selected_code

            with open("parameters/bands.yaml") as file:
                builder.bands.bands.pw.parameters = Dict(
                    dict=yaml.load(file, Loader=yaml.FullLoader)
                )
            builder.bands.bands.pw.metadata.options = self.options
            builder.bands.bands.pseudo_family = Str(self.pseudo_family_selector.value)
            builder.bands.bands.pw.code = self.code_group_pw.selected_code

        # PDOS section.
        if self.run_pdos.value:
            pass

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
