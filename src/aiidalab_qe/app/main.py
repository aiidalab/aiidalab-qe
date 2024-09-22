"""The main widget that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.configuration.model import ConfigurationModel
from aiidalab_qe.app.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.result.model import ResultsModel
from aiidalab_qe.app.structure import StructureSelectionStep
from aiidalab_qe.app.structure.model import StructureModel
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_qe.app.submission.model import SubmissionModel
from aiidalab_qe.common import QeAppWorkChainSelector
from aiidalab_widgets_base import WizardAppWidget


class App(ipw.VBox):
    """The main widget that combines all the application steps together."""

    def __init__(self, qe_auto_setup=True):
        # Initialize the models
        struct_model = StructureModel()
        config_model = ConfigurationModel()
        submit_model = SubmissionModel()
        results_model = ResultsModel()

        # Create the application steps
        self.structure_step = StructureSelectionStep(
            model=struct_model,
            auto_advance=True,
        )

        self.configure_step = ConfigureQeAppWorkChainStep(
            model=config_model,
            auto_advance=True,
        )

        self.submit_step = SubmitQeAppWorkChainStep(
            model=submit_model,
            auto_advance=True,
            qe_auto_setup=qe_auto_setup,
        )

        self.results_step = ViewQeAppWorkChainStatusAndResultsStep(model=results_model)

        # Link the models of the application steps
        ipw.dlink(
            (self.structure_step, "state"),
            (self.configure_step, "previous_step_state"),
        )
        ipw.dlink(
            (struct_model, "confirmed_structure"),
            (config_model, "input_structure"),
        )
        ipw.dlink(
            (struct_model, "confirmed_structure"),
            (submit_model, "input_structure"),
        )
        ipw.dlink(
            (self.configure_step, "state"),
            (self.submit_step, "previous_step_state"),
        )
        ipw.dlink(
            (config_model, "configuration_parameters"),
            (submit_model, "input_parameters"),
        )
        ipw.dlink(
            (submit_model, "process"),
            (results_model, "process"),
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
        self._wizard_app_widget.observe(
            self._observe_selected_index,
            "selected_index",
        )

        # Add process selection header
        self.work_chain_selector = QeAppWorkChainSelector(
            layout=ipw.Layout(width="auto")
        )
        self.work_chain_selector.observe(
            self._observe_process_selection,
            "value",
        )

        ipw.dlink(
            (self.submit_step, "process"),
            (self.work_chain_selector, "value"),
            transform=lambda node: None if node is None else node.pk,
        )

        super().__init__(
            children=[
                self.work_chain_selector,
                self._wizard_app_widget,
            ]
        )

        self._wizard_app_widget.selected_index = None

        # TODO disable other panels until required components are rendered

    @property
    def steps(self):
        return self._wizard_app_widget.steps

    def _observe_selected_index(self, change):
        """Check unsaved change in the step when leaving the step."""
        from aiidalab_widgets_base import WizardAppWidgetStep

        # no accordion tab is selected
        if change["new"] is None:
            return
        new_idx = change["new"]

        step = self.steps[new_idx][1]
        step.render()

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

    def _observe_process_selection(self, change):
        from aiida.orm import load_node
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
