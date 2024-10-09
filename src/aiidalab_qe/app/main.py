"""The main widget that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

import ipywidgets as ipw
import traitlets as tl
from IPython.display import Javascript, display

from aiida.orm import load_node
from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.structure import StructureSelectionStep
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_widgets_base import WizardAppWidget, WizardAppWidgetStep


class App(ipw.VBox):
    """The main widget that combines all the application steps together."""

    # The PK or UUID of the work chain node.
    process = tl.Union([tl.Unicode(), tl.Int()], allow_none=True)

    def __init__(self, qe_auto_setup=True):
        # Create the application steps
        self.structure_step = StructureSelectionStep(auto_advance=True)
        self.structure_step.observe(self._observe_structure_selection, "structure")
        self.configure_step = ConfigureQeAppWorkChainStep(auto_advance=True)
        self.submit_step = SubmitQeAppWorkChainStep(
            auto_advance=True,
            qe_auto_setup=qe_auto_setup,
        )
        self.results_step = ViewQeAppWorkChainStatusAndResultsStep()

        # Link the application steps
        ipw.dlink(
            (self.structure_step, "state"),
            (self.configure_step, "previous_step_state"),
        )
        ipw.dlink(
            (self.structure_step, "confirmed_structure"),
            (self.submit_step, "input_structure"),
        )
        ipw.dlink(
            (self.structure_step, "confirmed_structure"),
            (self.configure_step, "input_structure"),
        )
        ipw.dlink(
            (self.configure_step, "state"),
            (self.submit_step, "previous_step_state"),
        )
        ipw.dlink(
            (self.configure_step, "configuration_parameters"),
            (self.submit_step, "input_parameters"),
        )
        ipw.dlink(
            (self.submit_step, "process"),
            (self.results_step, "process"),
            transform=lambda node: node.uuid if node is not None else None,
        )

        # Add the application steps to the application
        self._wizard_app_widget = WizardAppWidget(
            steps=[
                ("Select structure", self.structure_step),
                ("Configure workflow", self.configure_step),
                ("Choose computational resources", self.submit_step),
                ("Status & Results", self.results_step),
            ]
        )
        # hide the header
        self._wizard_app_widget.children[0].layout.display = "none"
        self._wizard_app_widget.observe(self._observe_selected_index, "selected_index")

        # Add a button to start a new calculation
        self.new_work_chains_button = ipw.Button(
            description="New calculation",
            tooltip="Start a new calculation",
            button_style="success",
            icon="plus-circle",
        )

        def on_button_click(_):
            display(Javascript("window.open('./qe.ipynb', '_blank')"))

        self.new_work_chains_button.on_click(on_button_click)

        super().__init__(
            children=[
                self.new_work_chains_button,
                self._wizard_app_widget,
            ]
        )

    @property
    def steps(self):
        return self._wizard_app_widget.steps

    # Reset the confirmed_structure in case that a new structure is selected
    def _observe_structure_selection(self, change):
        with self.structure_step.hold_sync():
            if (
                self.structure_step.confirmed_structure is not None
                and self.structure_step.confirmed_structure != change["new"]
            ):
                self.structure_step.confirmed_structure = None

    def _observe_selected_index(self, change):
        """Check unsaved change in the step when leaving the step."""
        # no accordion tab is selected
        if not change["new"]:
            return
        new_idx = change["new"]
        # only when entering the submit step, check and udpate the blocker messages
        # steps[new_idx][0] is the title of the step
        if self.steps[new_idx][1] is not self.submit_step:
            return
        blockers = []
        # Loop over all steps before the submit step
        for title, step in self.steps[:new_idx]:
            # check if the step is saved
            if not step.is_saved():
                step.state = WizardAppWidgetStep.State.CONFIGURED
                blockers.append(
                    f"Unsaved changes in the <b>{title}</b> step. Please save the changes before submitting."
                )
        self.submit_step.external_submission_blockers = blockers

    @tl.observe("process")
    def _observe_process(self, change):
        from aiida.orm.utils.serialize import deserialize_unsafe

        if change["old"] == change["new"]:
            return
        pk = change["new"]
        if pk is None:
            self._wizard_app_widget.reset()
            self._wizard_app_widget.selected_index = 0
        else:
            process = load_node(pk)
            with self.structure_step.manager.hold_sync():
                with self.structure_step.hold_sync():
                    self._wizard_app_widget.selected_index = 3
                    self.structure_step.manager.viewer.structure = (
                        process.inputs.structure.get_ase()
                    )
                    self.structure_step.structure = process.inputs.structure
                    self.structure_step.confirm()
                    self.submit_step.process = process

            # set ui_parameters
            # print out error message if yaml format ui_parameters is not reachable
            ui_parameters = process.base.extras.get("ui_parameters", {})
            if ui_parameters and isinstance(ui_parameters, str):
                ui_parameters = deserialize_unsafe(ui_parameters)
                self.configure_step.set_configuration_parameters(ui_parameters)
                self.configure_step.confirm()
                self.submit_step.set_submission_parameters(ui_parameters)
                self.submit_step.state = self.submit_step.State.SUCCESS
