"""The main widget that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

import ipywidgets as ipw
import traitlets as tl

from aiida.orm import load_node
from aiida.orm.utils.serialize import deserialize_unsafe
from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.configuration.model import ConfigurationModel
from aiidalab_qe.app.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.result.model import ResultsModel
from aiidalab_qe.app.structure import StructureSelectionStep
from aiidalab_qe.app.structure.model import StructureModel
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_qe.app.submission.model import SubmissionModel
from aiidalab_qe.common import QeAppWorkChainSelector
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base import WizardAppWidget, WizardAppWidgetStep


class App(ipw.VBox):
    """The main widget that combines all the application steps together."""

    blockers = tl.List(tl.Unicode(), default_value=[])

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

        # Link wizard steps
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
        ipw.dlink(
            (self, "blockers"),
            (submit_model, "external_submission_blockers"),
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
            self._on_step_selection,
            "selected_index",
        )

        # Add process selection header
        self.work_chain_selector = QeAppWorkChainSelector(
            layout=ipw.Layout(width="auto")
        )
        self.work_chain_selector.observe(
            self._on_process_selection,
            "value",
        )

        self._process_loading_message = LoadingWidget(
            message="Loading process",
            layout=ipw.Layout(display="none"),
        )

        ipw.dlink(
            (submit_model, "process"),
            (self.work_chain_selector, "value"),
            transform=lambda node: None if node is None else node.pk,
        )

        super().__init__(
            children=[
                self.work_chain_selector,
                self._process_loading_message,
                self._wizard_app_widget,
            ]
        )

        self._wizard_app_widget.selected_index = None

        # TODO disable other panels until required components are rendered

    @property
    def steps(self):
        return self._wizard_app_widget.steps

    def _on_step_selection(self, change):
        if (step_index := change["new"]) is None:
            return

        self._render_step(step_index)

        if self.steps[step_index][1] is not self.submit_step:
            return

        self._update_blockers(step_index)

    def _on_process_selection(self, change):
        self._update_from_process(change["new"])

    def _render_step(self, step_index):
        step = self.steps[step_index][1]
        step.render()

    def _update_blockers(self, step_index):
        with self.hold_trait_notifications():
            self.blockers = []
            for title, step in self.steps[:step_index]:
                if not step.is_saved():
                    step.state = WizardAppWidgetStep.State.CONFIGURED
                    self.blockers.append(
                        f"Unsaved changes in the <b>{title}</b> step. Please save the changes before submitting."
                    )

    def _update_from_process(self, pk):
        if pk is None:
            self._wizard_app_widget.reset()
            self._wizard_app_widget.selected_index = 0
        else:
            self._show_process_loading_message()
            process = load_node(pk)
            self.structure_step.render()
            self.configure_step.render()
            self.submit_step.render()
            with self.structure_step.manager.hold_sync():
                with self.structure_step.hold_sync():
                    self._wizard_app_widget.selected_index = 3
                    self.structure_step.set_structure(process.inputs.structure)
                    self.structure_step.confirm()
                    self.submit_step.set_process(process)
            parameters = process.base.extras.get("ui_parameters", {})
            if parameters and isinstance(parameters, str):
                parameters = deserialize_unsafe(parameters)
                self.configure_step.set_configuration_parameters(parameters)
                self.configure_step.confirm()
                self.submit_step.set_submission_parameters(parameters)
                self.submit_step.state = WizardAppWidgetStep.State.SUCCESS
            self._hide_process_loading_message()

    def _show_process_loading_message(self):
        self._process_loading_message.layout.display = "flex"

    def _hide_process_loading_message(self):
        self._process_loading_message.layout.display = "none"
