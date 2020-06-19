"""Widgets for the submission of codes.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import threading
from time import sleep

from IPython.display import display, clear_output
import ipywidgets as ipw
import traitlets
from aiida.orm import load_node, ProcessNode
from aiida.engine import ProcessState

from aiidalab_widgets_base import CodeDropdown, viewer

from process import ProcessStatusWidget
from wizard import WizardApp, WizardAppStep


class CodeSubmitWidget(ipw.VBox, WizardAppStep):

    process = traitlets.Instance(ProcessNode, allow_none=True)
    disabled = traitlets.Bool()

    def _update_total_num_cpus(self, change):
        self.total_num_cpus.value = self.number_of_nodes.value * self.cpus_per_node.value

    def __init__(self, description=None, **kwargs):
        self.code_group = CodeDropdown(input_plugin='quantumespresso.pw', text="Select code")

        extra = {
            'style': {'description_width': '150px'},
            'layout': {'max_width': '200px'}
        }

        self.number_of_nodes = ipw.BoundedIntText(
            value=1, step=1, min=1,
            description="# nodes",
            disabled=False,
            **extra)
        self.cpus_per_node = ipw.BoundedIntText(
            value=1, step=1, min=1,
            description="# cpus per node",
            **extra)
        self.total_num_cpus = ipw.BoundedIntText(
            value=1, step=1, min=1,
            description="# total cpus",
            disabled=True,
            **extra)

        # Update the total # of CPUs int text:
        self.number_of_nodes.observe(self._update_total_num_cpus, 'value')
        self.cpus_per_node.observe(self._update_total_num_cpus, 'value')

        resource_selection_prompt = ipw.HTML(
            "Select the compute resources for this calculation.")
        resource_selection_help = ipw.HTML("""<div style="line-height:120%; padding-top:25px;">
            <p>There is no general rule of thumb on how to select the appropriate number of
            nodes and cores. In general:</p>
            <ul>
            <li>Increase the number of nodes if you run out of memory for larger structures.</li>
            <li>Increase the number of nodes and cores if you want to reduce the total runtime.</li>
            </ul>
            <p>However, specifying the optimal configuration of resources is a complex issue and
            simply increasing either cores or nodes may not have the desired effect.</p></div>""")
        self.resources = ipw.HBox(children=[
            ipw.VBox(
                children=[
                    resource_selection_prompt,
                    self.number_of_nodes,
                    self.cpus_per_node,
                    self.total_num_cpus,
                ],
                layout=ipw.Layout(min_width='310px'),
            ),
            resource_selection_help,
        ])

        # Setup pseudofamily potential selection group:
        link_url = "https://www.materialscloud.org/discover/sssp/table/precision"
        pseudo_family_prompt = ipw.HTML(
            f'Select the <a href="{link_url}" target="_blank">pseudopotential library</a> for the calculation.')
        pseudo_family_help = ipw.HTML("""
            <div style="line-height:120%;">If you are unsure what to choose, select 'SSSP efficiency', which for most
            calculations will produce sufficiently accurate results at comparatively small computational cost. If
            your calculation requires a higher accuracy, select 'SSSP accuracy', which will be computationally more
            expensive, but will produce even more accurate results.</div>""")
        self.pseudo_family_selection = ipw.ToggleButtons(
            options={
                'SSSP efficiency': 'SSSP_1.1_efficiency',
                'SSSP accuracy': 'SSSP_1.1_precision',
            },
        )
        self.pseudo_family_group = ipw.VBox(children=[
            pseudo_family_prompt,
            self.pseudo_family_selection,
            pseudo_family_help,
            ])

        # Clicking on the 'submit' button will trigger the execution of the
        # submit() method.
        self.submit_button = ipw.Button(
            description='Submit',
            tooltip="Submit the calculation with the selected parameters.",
            icon='play',
            button_style='success',
            layout=ipw.Layout(width='auto', flex="1 1 auto"),
            disabled=True)
        self.submit_button.on_click(self._on_submit_button_clicked)

        # The 'skip' button is only shown when the skip() method is implemented.
        self.skip_button = ipw.Button(
            description='Skip',
            icon='fast-forward',
            button_style='info',
            layout=ipw.Layout(width='auto', flex="1 1 auto"),
            disabled=True)
        if self.skip:  # skip() method is implemented
            # connect with skip_button
            self.skip_button.on_click(self.skip)  # connect with skip_button
        else:  # skip() not implemented
            # hide the button
            self.skip_button.layout.visibility = 'hidden'

        # Place all buttons at the footer of the widget.
        self.buttons = ipw.HBox(children=[self.submit_button, self.skip_button])

        self.config_tabs = ipw.Tab(
            children=[self.code_group, self.pseudo_family_group, self.resources],
            layout=ipw.Layout(height='200px'),
        )
        self.config_tabs.set_title(0, 'Code')
        self.config_tabs.set_title(1, 'Pseudopotential')
        self.config_tabs.set_title(2, 'Compute resources')

        self.process_status = ProcessStatusWidget()
        ipw.dlink((self, 'process'), (self.process_status, 'process'))

        self.outputs_keys = ipw.Dropdown()
        self.outputs_keys.observe(self._refresh_outputs_view, names=['options', 'value'])
        self.output_control = ipw.HBox(children=[self.outputs_keys])

        self.output_area = ipw.Output(
            layout={
                'width': 'auto',
                'height': 'auto',
                'border': '1px solid black'})
        self.results_view = ipw.VBox(children=[self.output_control, self.output_area])

        self.accordion = ipw.Accordion(
            children=[self.config_tabs, self.process_status, self.results_view])
        self.accordion.set_title(0, 'Config')
        self.accordion.set_title(1, 'Status')
        self.accordion.set_title(2, 'Results (0)')

        self.callbacks = list()

        ipw.dlink((self, 'disabled'), (self.skip_button, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.submit_button, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.code_group.dropdown, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.number_of_nodes, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.cpus_per_node, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.pseudo_family_selection, 'disabled'))

        # Initialize widget disabled status based on step state.
        self.disabled = self.state != WizardApp.State.READY

        super().__init__(children=[
            ipw.Label() if description is None else description,
            self.accordion,
            self.buttons],
            **kwargs)

    def _update_state(self):
        "Update state based on the process state."
        if self.process is None:
            self.state = WizardApp.State.INIT
        else:
            process_state = load_node(self.process.id).process_state
            if process_state in (ProcessState.CREATED, ProcessState.RUNNING, ProcessState.WAITING):
                self.state = WizardApp.State.ACTIVE
            elif process_state in (ProcessState.EXCEPTED, ProcessState.KILLED):
                self.state = WizardApp.State.FAIL
            elif process_state is ProcessState.FINISHED:
                self.state = WizardApp.State.SUCCESS

    @traitlets.observe('state')
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.disabled = change['new'] != WizardApp.State.READY

            if change['new'] == WizardApp.State.ACTIVE:
                self.accordion.selected_index = 1

    @property
    def options(self):
        return {
            'max_wallclock_seconds': 3600*2,
            'resources': {
                'num_machines': self.number_of_nodes.value,
                'num_mpiprocs_per_machine': self.cpus_per_node.value}}

    def _refresh_outputs_keys(self):
        if self.process is None:
            with self.hold_sync():
                self.outputs_keys.options = ["[No outputs]"]
                self.outputs_keys.disabled = True

        else:
            with self.hold_sync():
                process_node = load_node(self.process.id)
                self.outputs_keys.options = ["[Select output]"] + \
                    [str(o) for o in process_node.outputs]
            return process_node

    def _refresh_outputs_view(self, change=None):
        self.accordion.set_title(2, f"Results ({len(self.outputs_keys.options)-1})")

        if change is None or change['name'] == 'options':
            selection_key = self.outputs_keys.value
        else:
            selection_key = change['new']

        with self.output_area:
            # Clear first to ensure that we are not showing the wrong thing.
            clear_output()

            if self.outputs_keys.index and selection_key is not None:
                # Load data
                process_node = load_node(self.process.id)
                data = process_node.outputs[selection_key]
                output_viewer_widget = viewer(data)

                display(output_viewer_widget)

    def _monitor_process(self):
        assert self.process is not None
        process_node = load_node(self.process.id)

        while not process_node.is_sealed:
            self.process_status.update()
            for callback in self.callbacks:
                callback(self)
            sleep(0.1)

        with self.hold_trait_notifications():
            self.process_status.update()
            self._refresh_outputs_keys()

        return process_node

    @traitlets.observe('process')
    def _observe_process(self, change):
        process = change['new']
        with self.hold_trait_notifications():
            if process is None:
                self._refresh_outputs_keys()
            else:
                process_node = load_node(process.id)
                self.process_status.process = process_node
                if process_node.is_sealed:
                    self._refresh_outputs_keys()
                else:
                    self.state = WizardApp.State.ACTIVE
                    monitor_thread = threading.Thread(target=self._monitor_process)
                    monitor_thread.start()
                return process_node

    skip = False

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.state = WizardApp.State.ACTIVE
        self.submit()

    def submit(self, _):
        raise NotImplementedError()
