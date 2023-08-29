# -*- coding: utf-8 -*-
"""The main widget that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

import ipywidgets as ipw
from aiida.orm import load_node
from aiidalab_widgets_base import (
    BasicCellEditor,
    BasicStructureEditor,
    OptimadeQueryWidget,
    StructureBrowserWidget,
    StructureExamplesWidget,
    StructureManagerWidget,
    StructureUploadWidget,
    WizardAppWidget,
    WizardAppWidgetStep,
)

from aiidalab_qe.app.common import AddingTagsEditor, QeAppWorkChainSelector
from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.structure import Examples, StructureSelectionStep
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep

OptimadeQueryWidget.title = "OPTIMADE"  # monkeypatch


class App(ipw.VBox):
    """The main widget that combines all the application steps together."""

    def __init__(self, qe_auto_setup=True):
        # Create the application steps
        self.structure_manager_widget = StructureManagerWidget(
            importers=[
                StructureUploadWidget(title="Upload file"),
                OptimadeQueryWidget(embedded=False),
                StructureBrowserWidget(title="AiiDA database"),
                StructureExamplesWidget(title="From Examples", examples=Examples),
            ],
            editors=[
                BasicCellEditor(title="Edit cell"),
                BasicStructureEditor(title="Edit structure"),
                AddingTagsEditor(title="Edit tags"),
            ],
            node_class="StructureData",
            storable=False,
            configuration_tabs=["Cell", "Selection", "Appearance", "Download"],
        )
        self.structure_selection_step = StructureSelectionStep(
            manager=self.structure_manager_widget, auto_advance=True
        )
        self.structure_selection_step.observe(
            self._observe_structure_selection, "structure"
        )

        self.configure_qe_app_work_chain_step = ConfigureQeAppWorkChainStep(
            auto_advance=True
        )
        self.submit_qe_app_work_chain_step = SubmitQeAppWorkChainStep(
            auto_advance=True, qe_auto_setup=qe_auto_setup
        )
        view_qe_app_work_chain_status_and_results_step = (
            ViewQeAppWorkChainStatusAndResultsStep()
        )

        # Link the application steps
        ipw.dlink(
            (self.structure_selection_step, "state"),
            (self.configure_qe_app_work_chain_step, "previous_step_state"),
        )
        ipw.dlink(
            (self.structure_selection_step, "confirmed_structure"),
            (self.submit_qe_app_work_chain_step, "input_structure"),
        )
        ipw.dlink(
            (self.structure_selection_step, "confirmed_structure"),
            (self.configure_qe_app_work_chain_step, "input_structure"),
        )
        ipw.dlink(
            (self.configure_qe_app_work_chain_step, "state"),
            (self.submit_qe_app_work_chain_step, "previous_step_state"),
        )
        ipw.dlink(
            (self.configure_qe_app_work_chain_step, "workchain_settings"),
            (self.submit_qe_app_work_chain_step, "workchain_settings"),
        )
        ipw.dlink(
            (self.configure_qe_app_work_chain_step, "advanced_settings"),
            (self.submit_qe_app_work_chain_step, "advanced_settings"),
        )
        ipw.dlink(
            (self.configure_qe_app_work_chain_step, "pseudo_family_selector"),
            (self.submit_qe_app_work_chain_step, "pseudo_family_selector"),
        )
        ipw.dlink(
            (self.configure_qe_app_work_chain_step, "pseudo_setter"),
            (self.submit_qe_app_work_chain_step, "pseudo_setter"),
        )

        ipw.dlink(
            (self.submit_qe_app_work_chain_step, "process"),
            (view_qe_app_work_chain_status_and_results_step, "process"),
            transform=lambda node: node.uuid if node is not None else None,
        )

        # Add the application steps to the application
        self._wizard_app_widget = WizardAppWidget(
            steps=[
                ("Select structure", self.structure_selection_step),
                ("Configure workflow", self.configure_qe_app_work_chain_step),
                ("Choose computational resources", self.submit_qe_app_work_chain_step),
                ("Status & Results", view_qe_app_work_chain_status_and_results_step),
            ]
        )

        # Add process selection header
        work_chain_selector = QeAppWorkChainSelector(layout=ipw.Layout(width="auto"))
        work_chain_selector.observe(self._observe_process_selection, "value")

        ipw.dlink(
            (self.submit_qe_app_work_chain_step, "process"),
            (work_chain_selector, "value"),
            transform=lambda node: None if node is None else node.pk,
        )

        super().__init__(
            children=[
                work_chain_selector,
                self._wizard_app_widget,
            ]
        )

    # Reset all subsequent steps in case that a new structure is selected
    def _observe_structure_selection(self, change):
        with self.structure_selection_step.hold_sync():
            if (
                self.structure_selection_step.confirmed_structure is not None
                and self.structure_selection_step.confirmed_structure != change["new"]
            ):
                self._wizard_app_widget.reset()

    def _observe_process_selection(self, change):
        if change["old"] == change["new"]:
            return
        pk = change["new"]
        if pk is None:
            self._wizard_app_widget.reset()
            self._wizard_app_widget.selected_index = 0
        else:
            process = load_node(pk)
            with self.structure_manager_widget.hold_sync():
                with self.structure_selection_step.hold_sync():
                    self._wizard_app_widget.selected_index = 3
                    self.structure_manager_widget.input_structure = (
                        process.inputs.structure
                    )
                    self.structure_selection_step.structure = process.inputs.structure
                    self.structure_selection_step.confirmed_structure = (
                        process.inputs.structure
                    )
                    self.configure_qe_app_work_chain_step.state = (
                        WizardAppWidgetStep.State.SUCCESS
                    )
                    self.submit_qe_app_work_chain_step.process = process
