# -*- coding: utf-8 -*-
"""The main widget that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

import copy
import typing as t

import ipywidgets as ipw
from aiida.orm import load_node
from aiidalab_widgets_base import WizardAppWidget, WizardAppWidgetStep

from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.structure import StructureSelectionStep
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_qe.common import QeAppWorkChainSelector


class App(ipw.VBox):
    """The main widget that combines all the application steps together."""

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
        self._wizard_app_widget.observe(self._observe_selected_index, "selected_index")

        # Add process selection header
        self.work_chain_selector = QeAppWorkChainSelector(
            layout=ipw.Layout(width="auto")
        )
        self.work_chain_selector.observe(self._observe_process_selection, "value")

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

    @property
    def steps(self) -> t.List[t.Tuple[str, WizardAppWidgetStep]]:
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
        if not change["new"]:
            # not accordion tab selected
            return

        new_idx = change["new"]
        # if entering the submit step, udpate the blocker messages
        _, step = self.steps[new_idx]

        if step is not self.submit_step:
            return

        # Loop over all steps before the submit step and check if they are ready.
        for title, preceding_step in self.steps[:new_idx]:
            # check if the step is saved
            blockers = copy.deepcopy(self.submit_step.external_submission_blockers)
            if not preceding_step.is_saved():
                preceding_step.state = WizardAppWidgetStep.State.CONFIGURED
                blockers.add(
                    f"Unsaved changes in the <b>{title}</b> step. Please save the changes before submitting."
                )
            else:
                # remove the blocker message if the step has all changes confirmed
                blockers = set()

            self.submit_step.external_submission_blockers = blockers

    def _observe_process_selection(self, change):
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
                    self.structure_step.confirmed_structure = process.inputs.structure
                    self.configure_step.state = WizardAppWidgetStep.State.SUCCESS
                    self.submit_step.process = process
            # set ui_parameters
            # print out error message if yaml format ui_parameters is not reachable
            ui_parameters = process.base.extras.get("ui_parameters", {})
            if ui_parameters and isinstance(ui_parameters, str):
                ui_parameters = deserialize_unsafe(ui_parameters)
                self.configure_step.set_configuration_parameters(ui_parameters)
                self.configure_step.state = self.configure_step.State.SUCCESS
                self.submit_step.set_submission_parameters(ui_parameters)
                self.submit_step.state = self.submit_step.State.SUCCESS
