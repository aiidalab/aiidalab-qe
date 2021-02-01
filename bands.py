"""Widgets for the submission of bands work chains.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import ipywidgets as ipw
import traitlets
from aiida.engine import ProcessState
from aiida.engine import submit
from aiidalab_widgets_base import CodeDropdown, viewer
from aiida.orm import load_node, ProcessNode
from aiida.orm import StructureData, Float, Str, BandsData
from aiida.plugins import WorkflowFactory
from IPython.display import display, clear_output

from codes import ResourceSelectionWidget
from process import ProcessMonitor
from process import ProcessStatusWidget
from pseudos import PseudoFamilySelector
from util import load_default_parameters
from wizard import WizardApp, WizardAppStep


WARNING_ICON = "\u26A0"


class ComputeBandsSubmitWidget(ipw.VBox, WizardAppStep):
    """Step for submission of a bands workchain."""

    input_structure = traitlets.Instance(StructureData, allow_none=True)
    band_structure = traitlets.Instance(BandsData, allow_none=True)

    process = traitlets.Instance(ProcessNode, allow_none=True)
    disabled = traitlets.Bool()

    def __init__(self, description=None, **kwargs):
        self.code_group = CodeDropdown(
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
        self.code_group.observe(lambda _: self._update_state(), ["selected_code"])

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

        # The 'skip' button is only shown when the skip() method is implemented.
        self.skip_button = ipw.Button(
            description="Skip",
            icon="fast-forward",
            button_style="info",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )
        if self.skip:  # skip() method is implemented
            # connect with skip_button
            self.skip_button.on_click(self.skip)  # connect with skip_button
        else:  # skip() not implemented
            # hide the button
            self.skip_button.layout.visibility = "hidden"

        # Place all buttons at the footer of the widget.
        self.buttons = ipw.HBox(children=[self.submit_button, self.skip_button])

        self.config_tabs = ipw.Tab(
            children=[self.code_group, self.pseudo_family_selector, self.resources],
            layout=ipw.Layout(height="200px"),
        )
        self.config_tabs.set_title(0, "Code")
        # second tab initialized below
        self.config_tabs.set_title(2, "Compute resources")

        # Show warning in cofig title when pseudos are not installed:
        def _observe_sssp_installed(change):
            self.config_tabs.set_title(
                1, "Pseudopotential" + ("" if change["new"] else f" {WARNING_ICON}")
            )
            self._observe_state(change=dict(new=self.state))  # trigger refresh

        self.pseudo_family_selector.observe(_observe_sssp_installed, "installed")
        _observe_sssp_installed(
            change=dict(new=self.pseudo_family_selector.installed)
        )  # init

        self.process_status = ProcessStatusWidget()
        ipw.dlink((self, "process"), (self.process_status, "process"))

        self.outputs_keys = ipw.Dropdown()
        self.outputs_keys.observe(
            self._refresh_outputs_view, names=["options", "value"]
        )
        self.output_area = ipw.Output(
            layout={"width": "auto", "height": "auto", "border": "1px solid black"}
        )
        self.results_view = ipw.VBox(children=[self.outputs_keys, self.output_area])

        self.accordion = ipw.Accordion(
            children=[self.config_tabs, self.process_status, self.results_view]
        )
        self.accordion.set_title(0, "Config")
        self.accordion.set_title(1, "Status")
        self.accordion.set_title(2, "Results (0)")

        ipw.dlink((self, "disabled"), (self.skip_button, "disabled"))
        ipw.dlink((self, "disabled"), (self.code_group.dropdown, "disabled"))
        ipw.dlink((self, "disabled"), (self.resources.number_of_nodes, "disabled"))
        ipw.dlink((self, "disabled"), (self.resources.cpus_per_node, "disabled"))
        ipw.dlink((self, "disabled"), (self.pseudo_family_selector, "disabled"))

        # Initialize widget disabled status based on step state.
        self.disabled = self.state != WizardApp.State.READY

        # Setup process monitor
        self.process_monitor = ProcessMonitor(
            callbacks=[
                (lambda _: self.process_status.update(), 1),
                (self._refresh_outputs_keys, 1000),
            ]
        )
        ipw.dlink((self, "process"), (self.process_monitor, "process"))

        super().__init__(
            children=[
                ipw.Label(
                    'Specify the parameters and options for the calculation and then click on "Submit".'
                ),
                self.accordion,
                self.buttons,
            ],
            **kwargs,
        )

    def _update_state(self):
        if self.process is None:
            if self.input_structure is None:
                if self.band_structure is None:
                    self.state = WizardApp.State.INIT
                else:
                    self.state = WizardApp.State.FAIL
            else:
                if self.band_structure is None:
                    if self.code_group.selected_code is None:
                        self.state = WizardApp.State.READY
                    else:
                        self.state = WizardApp.State.CONFIGURED
                else:
                    self.state = WizardApp.State.SUCCESS
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

    @traitlets.observe("state")
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.disabled = (
                change["new"] not in (WizardApp.State.READY, WizardApp.State.CONFIGURED)
                or not self.pseudo_family_selector.installed
            )
            self.submit_button.disabled = change["new"] != WizardApp.State.CONFIGURED

            if change["new"] == WizardApp.State.ACTIVE:
                self.accordion.selected_index = 1

    @traitlets.observe("input_structure")
    def _observe_input_structure(self, change):
        self._update_state()

    @traitlets.observe("band_structure")
    def _observe_band_structure(self, change):
        self._update_state()

    @traitlets.validate("process")
    def _validate_process(self, proposal):
        process_node = proposal["value"]
        builder = process_node.get_builder_restart()
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

        return proposal["value"]

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is None:
                self.band_structure = None
            else:
                # Restore the code parameters
                builder = process_node.get_builder_restart()
                self.pseudo_family_selector.value = builder.scf.pseudo_family.value
                try:
                    self.code_group.selected_code = builder.scf.pw["code"].full_label
                except traitlets.TraitError:
                    # It's not considered a critical issue if the exact code cannot
                    # be retrieved, however we set it to "None" in the interface to
                    # indicate that to the user.
                    self.code_group.selected_code = None
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

    def _refresh_outputs_keys(self, process_id):
        with self.hold_trait_notifications():
            if process_id is None:
                self.outputs_keys.options = ["[No outputs]"]
                self.outputs_keys.disabled = True
                return None

            process_node = load_node(process_id)
            self.outputs_keys.options = ["[Select output]"] + [
                str(o) for o in process_node.outputs
            ]
            if "band_structure" in process_node.outputs:
                self.band_structure = process_node.outputs.band_structure
            self.outputs_keys.disabled = False
            return process_node

    def _refresh_outputs_view(self, change=None):
        self.accordion.set_title(2, f"Results ({len(self.outputs_keys.options)-1})")

        if change is None or change["name"] == "options":
            selection_key = self.outputs_keys.value
        else:
            selection_key = change["new"]

        with self.output_area:
            # Clear first to ensure that we are not showing the wrong thing.
            clear_output()

            if self.outputs_keys.index and selection_key is not None:
                # Load data
                process_node = load_node(self.process.id)
                data = process_node.outputs[selection_key]
                output_viewer_widget = viewer(data)

                display(output_viewer_widget)

    skip = False

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.state = WizardApp.State.ACTIVE
        self.submit()

    def submit(self, _=None):
        assert self.input_structure is not None

        builder = WorkflowFactory("quantumespresso.pw.bands").get_builder()

        builder.scf.pw.code = self.code_group.selected_code
        builder.scf.pw.parameters = load_default_parameters()
        builder.scf.pw.metadata.options = self.options
        builder.scf.kpoints_distance = Float(0.8)
        builder.scf.pseudo_family = Str(self.pseudo_family_selector.value)

        builder.bands.pw.code = self.code_group.selected_code
        builder.bands.pw.parameters = load_default_parameters()
        builder.bands.pw.metadata.options = self.options
        builder.bands.pseudo_family = Str(self.pseudo_family_selector.value)

        builder.structure = self.input_structure

        self.process = submit(builder)

    def reset(self):
        self.process = None
