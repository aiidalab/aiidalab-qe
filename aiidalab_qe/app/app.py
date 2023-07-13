import ipywidgets as ipw
from aiidalab_widgets_base import WizardAppWidget

from aiidalab_qe.app.configure.configure import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.process import QeAppWorkChainSelector
from aiidalab_qe.app.result.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.structures import StructureSelectionStep, structure_manager_widget
from aiidalab_qe.app.submit.submit import SubmitQeAppWorkChainStep


class QEApp:
    def __init__(self, qe_auto_setup=True) -> None:
        # To keep the Qeapp simple and clear, I moved the
        # structure_manager_widget to structure.py
        # Create the application steps
        self.structure_step = StructureSelectionStep(
            parent=self, manager=structure_manager_widget, auto_advance=True
        )
        self.configure_step = ConfigureQeAppWorkChainStep(
            parent=self, auto_advance=True
        )
        self.submit_step = SubmitQeAppWorkChainStep(
            parent=self, auto_advance=True, qe_auto_setup=qe_auto_setup
        )
        self.results_step = ViewQeAppWorkChainStatusAndResultsStep(parent=self)
        # Add the application steps to the application
        self.steps = WizardAppWidget(
            steps=[
                ("Select structure", self.structure_step),
                ("Configure workflow", self.configure_step),
                ("Choose computational resources", self.submit_step),
                ("Status & Results", self.results_step),
            ]
        )
        # Link the application steps
        ipw.dlink(
            (self.structure_step, "state"),
            (self.configure_step, "previous_step_state"),
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
            (self.structure_step, "confirmed_structure"),
            (self.submit_step, "input_structure"),
        )
        ipw.dlink(
            (self.submit_step, "process"),
            (self.results_step, "process"),
            transform=lambda node: node.uuid if node is not None else None,
        )

        # Reset all subsequent steps in case that a new structure is selected
        def _observe_structure_selection(change):
            with self.structure_step.hold_sync():
                if (
                    self.structure_step.confirmed_structure is not None
                    and self.structure_step.confirmed_structure != change["new"]
                ):
                    self.steps.reset()

        self.structure_step.observe(_observe_structure_selection, "structure")

        # Add process selection header
        self.work_chain_selector = QeAppWorkChainSelector(
            layout=ipw.Layout(width="auto")
        )

        def _observe_process_selection(change):
            from aiida.orm import load_node
            from aiidalab_widgets_base import WizardAppWidgetStep

            if change["old"] == change["new"]:
                return
            pk = change["new"]
            if pk is None:
                self.steps.reset()
                self.steps.selected_index = 0
                self.configure_step.reset()
                self.submit_step.reset()
            else:
                process = load_node(pk)
                with structure_manager_widget.hold_sync():
                    with self.structure_step.hold_sync():
                        self.steps.selected_index = 3
                        structure_manager_widget.input_structure = (
                            process.inputs.structure
                        )
                        self.structure_step.structure = process.inputs.structure
                        self.structure_step.confirmed_structure = (
                            process.inputs.structure
                        )
                        self.configure_step.state = WizardAppWidgetStep.State.SUCCESS
                        self.submit_step.process = process
                # set ui_parameters
                ui_parameters = process.base.extras.get("ui_parameters", {})
                self.configure_step.set_input_parameters(ui_parameters)
                self.submit_step.set_resources(ui_parameters["resources"])

        self.work_chain_selector.observe(_observe_process_selection, "value")
        ipw.dlink(
            (self.submit_step, "process"),
            (self.work_chain_selector, "value"),
            transform=lambda node: None if node is None else node.pk,
        )
        self.work_chain = ipw.VBox(children=[self.work_chain_selector, self.steps])
