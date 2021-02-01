"""Widgets for the submission of codes.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import itertools
import threading

from IPython.display import display, clear_output
import ipywidgets as ipw
import traitlets
from aiida.orm import load_node, ProcessNode
from aiida.engine import ProcessState

from aiidalab_widgets_base import CodeDropdown, viewer

from process import ProcessStatusWidget
from wizard import WizardApp, WizardAppStep
from pseudos import PseudoFamilySelector


WARNING_ICON = '\u26A0'


class ResourceSelectionWidget(ipw.HBox):
    """Widget for the selection of compute (CPU) resources."""

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


    def __init__(self, **kwargs):
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

        super().__init__(children=[
            ipw.VBox(
                children=[
                    self.resource_selection_prompt,
                    self.number_of_nodes,
                    self.cpus_per_node,
                    self.total_num_cpus,
                ],
                layout=ipw.Layout(min_width='310px'),
            ),
            self.resource_selection_help,
        ])

    def _update_total_num_cpus(self, change):
        self.total_num_cpus.value = self.number_of_nodes.value * self.cpus_per_node.value


class CodeSubmitWidget(ipw.VBox, WizardAppStep):
    """A step that concludes with the submission of a code."""

    process = traitlets.Instance(ProcessNode, allow_none=True)
    disabled = traitlets.Bool()

    def __init__(self, description=None, **kwargs):
        setup_code_params = {
            "computer": "localhost",
            "description":  "pw.x in AiiDAlab container.",
            "label": "pw",
            "input_plugin": "quantumespresso.pw",
            'remote_abs_path': '/usr/bin/pw.x',
        }
        self.code_group = CodeDropdown(input_plugin='quantumespresso.pw', text="Select code", setup_code_params=setup_code_params)
        self.code_group.observe(lambda _: self._update_state(), ['selected_code'])

        # Setup pseudo potential family selection
        self.pseudo_family_selector = PseudoFamilySelector()

        # Setup the compute resources tab
        self.resources = ResourceSelectionWidget()

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
            children=[self.code_group, self.pseudo_family_selector, self.resources],
            layout=ipw.Layout(height='200px'),
        )
        self.config_tabs.set_title(0, 'Code')
        # second tab initialized below
        self.config_tabs.set_title(2, 'Compute resources')

        # Show warning in cofig title when pseudos are not installed:
        def _observe_sssp_installed(change):
            self.config_tabs.set_title(1, 'Pseudopotential' + ('' if change['new'] else f' {WARNING_ICON}'))
            self._observe_state(change=dict(new=self.state))  # trigger refresh

        self.pseudo_family_selector.observe(_observe_sssp_installed, 'installed')
        _observe_sssp_installed(change=dict(new=self.pseudo_family_selector.installed))  # init

        self.process_status = ProcessStatusWidget()
        ipw.dlink((self, 'process'), (self.process_status, 'process'))

        self.outputs_keys = ipw.Dropdown()
        self.outputs_keys.observe(self._refresh_outputs_view, names=['options', 'value'])
        self.output_area = ipw.Output(
            layout={
                'width': 'auto',
                'height': 'auto',
                'border': '1px solid black'})
        self.results_view = ipw.VBox(children=[self.outputs_keys, self.output_area])

        self.accordion = ipw.Accordion(
            children=[self.config_tabs, self.process_status, self.results_view])
        self.accordion.set_title(0, 'Config')
        self.accordion.set_title(1, 'Status')
        self.accordion.set_title(2, 'Results (0)')

        self.callbacks = list()

        ipw.dlink((self, 'disabled'), (self.skip_button, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.code_group.dropdown, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.resources.number_of_nodes, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.resources.cpus_per_node, 'disabled'))
        ipw.dlink((self, 'disabled'), (self.pseudo_family_selector, 'disabled'))

        # Initialize widget disabled status based on step state.
        self.disabled = self.state != WizardApp.State.READY

        # Set up process monitoring.
        self._monitor_thread = None
        self._monitor_thread_stop = threading.Event()
        self._monitor_thread_lock = threading.Lock()

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
            process_state = self.process.process_state
            if process_state in (ProcessState.CREATED, ProcessState.RUNNING, ProcessState.WAITING):
                self.state = WizardApp.State.ACTIVE
            elif process_state in (ProcessState.EXCEPTED, ProcessState.KILLED):
                self.state = WizardApp.State.FAIL
            elif process_state is ProcessState.FINISHED:
                self.state = WizardApp.State.SUCCESS

    @traitlets.observe('state')
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.disabled = change['new'] not in (WizardApp.State.READY, WizardApp.State.CONFIGURED) \
                    or not self.pseudo_family_selector.installed
            self.submit_button.disabled = change['new'] != WizardApp.State.CONFIGURED

            if change['new'] == WizardApp.State.ACTIVE:
                self.accordion.selected_index = 1

    @property
    def options(self):
        return {
            'max_wallclock_seconds': 3600*2,
            'resources': {
                'num_machines': self.resources.number_of_nodes.value,
                'num_mpiprocs_per_machine': self.resources.cpus_per_node.value}}

    def _refresh_outputs_keys(self, process_id):
        with self.hold_trait_notifications():
            if process_id is None:
                self.outputs_keys.options = ["[No outputs]"]
                self.outputs_keys.disabled = True
                return None

            process_node = load_node(process_id)
            self.outputs_keys.options = ["[Select output]"] + \
                [str(o) for o in process_node.outputs]
            self.outputs_keys.disabled = False
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

    def _monitor_process(self, process_id, timeout=.1):
        self._monitor_thread_stop.wait(timeout=10 * timeout)  # brief delay to increase app stability

        process = None if process_id is None else load_node(process_id)

        iteration = itertools.count()
        while not (process is None or process.is_sealed):
            if next(iteration) % round(100 / timeout) == 0:
                self._refresh_outputs_keys(process_id)

            self.process_status.update()
            for callback in self.callbacks:
                callback(self)

            if self._monitor_thread_stop.wait(timeout=timeout):
                break

        # Final update:
        self.process_status.update()
        self._refresh_outputs_keys(process_id)

    @traitlets.observe('process')
    def _observe_process(self, change):
        self._update_state()
        process = change['new']
        if process is None or process.id != getattr(change['old'], 'id', None):
            with self.hold_trait_notifications():
                with self._monitor_thread_lock:
                    # stop thread
                    if self._monitor_thread is not None:
                        self._monitor_thread_stop.set()
                        self._monitor_thread.join()

                    # reset output
                    self._refresh_outputs_keys(None)

                    # start monitor thread
                    self._monitor_thread_stop.clear()
                    process_id = getattr(process, 'id', None)
                    self._monitor_thread = threading.Thread(target=self._monitor_process, args=(process_id, ))
                    self._monitor_thread.start()

    skip = False

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.state = WizardApp.State.ACTIVE
        self.submit()

    def submit(self, _):
        raise NotImplementedError()
