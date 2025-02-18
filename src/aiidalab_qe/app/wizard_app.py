import ipywidgets as ipw
import traitlets as tl
from IPython.display import Javascript, display

from aiida.orm import load_node
from aiida.orm.utils.serialize import deserialize_unsafe
from aiidalab_qe.app.configuration import ConfigureQeAppWorkChainStep
from aiidalab_qe.app.configuration.model import ConfigurationStepModel
from aiidalab_qe.app.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.result.model import ResultsStepModel
from aiidalab_qe.app.structure import StructureSelectionStep
from aiidalab_qe.app.structure.model import StructureStepModel
from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep
from aiidalab_qe.app.submission.model import SubmissionStepModel
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_qe.common.wizard import QeWizardStep
from aiidalab_widgets_base import WizardAppWidget


class WizardApp(ipw.VBox):
    """The main widget that combines all the application steps together."""

    # The PK or UUID of the work chain node.
    process = tl.Union([tl.Unicode(), tl.Int()], allow_none=True)

    def __init__(self, auto_setup=True, **kwargs):
        # Initialize the models
        self.structure_model = StructureStepModel()
        self.configure_model = ConfigurationStepModel()
        self.submit_model = SubmissionStepModel()
        self.results_model = ResultsStepModel()

        log_widget = kwargs.pop("log_widget", None)

        # Create the application steps
        self.structure_step = StructureSelectionStep(
            model=self.structure_model,
            auto_advance=True,
            auto_setup=auto_setup,
        )
        self.configure_step = ConfigureQeAppWorkChainStep(
            model=self.configure_model,
            auto_advance=True,
        )
        self.submit_step = SubmitQeAppWorkChainStep(
            model=self.submit_model,
            auto_advance=True,
            auto_setup=auto_setup,
        )
        self.results_step = ViewQeAppWorkChainStatusAndResultsStep(
            model=self.results_model,
            log_widget=log_widget,
        )

        # Wizard step observations
        ipw.dlink(
            (self.structure_step, "state"),
            (self.configure_step, "previous_step_state"),
        )
        self.structure_model.observe(
            self._on_structure_confirmation_change,
            "confirmed",
        )
        ipw.dlink(
            (self.configure_step, "state"),
            (self.submit_step, "previous_step_state"),
        )
        self.configure_model.observe(
            self._on_configuration_confirmation_change,
            "confirmed",
        )
        ipw.dlink(
            (self.submit_step, "state"),
            (self.results_step, "previous_step_state"),
        )
        self.submit_model.observe(
            self._on_submission,
            "confirmed",
        )

        # Add the application steps to the application
        self._wizard_app_widget = WizardAppWidget(
            steps=[
                ("Select structure", self.structure_step),
                ("Configure workflow", self.configure_step),
                ("Choose computational resources", self.submit_step),
                ("Status & results", self.results_step),
            ]
        )
        self._wizard_app_widget.observe(
            self._on_step_change,
            "selected_index",
        )

        # Hide the header
        self._wizard_app_widget.children[0].layout.display = "none"  # type: ignore

        self._process_loading_message = LoadingWidget(
            message="Loading process",
            layout={"display": "none"},
        )

        super().__init__(
            children=[
                InAppGuide(identifier="guide-header"),
                self._process_loading_message,
                self._wizard_app_widget,
                InAppGuide(identifier="post-guide"),
            ],
            **kwargs,
        )

        self._wizard_app_widget.selected_index = None

        self.structure_step.state = QeWizardStep.State.READY

    @property
    def steps(self):
        return self._wizard_app_widget.steps

    @tl.observe("process")
    def _on_process_change(self, change):
        self._update_from_process(change["new"])

    def _on_new_workchain_button_click(self, _):
        display(Javascript("window.open('./qe.ipynb', '_blank')"))

    def _on_step_change(self, change):
        if (step_index := change["new"]) is not None:
            self._render_step(step_index)

    def _on_structure_confirmation_change(self, _):
        self._update_configuration_step()

    def _on_configuration_confirmation_change(self, _):
        self._update_submission_step()

    def _on_submission(self, _):
        self._update_results_step()
        self._lock_app()

    def _render_step(self, step_index):
        step: QeWizardStep = self.steps[step_index][1]
        step.render()

    def _update_configuration_step(self):
        if self.structure_model.confirmed:
            self.configure_model.input_structure = self.structure_model.input_structure
        else:
            self.configure_model.input_structure = None

    def _update_submission_step(self):
        if self.configure_model.confirmed:
            self.submit_model.input_structure = self.structure_model.input_structure
            self.submit_model.input_parameters = self.configure_model.get_model_state()
        else:
            self.submit_model.input_structure = None
            self.submit_model.input_parameters = {}

    def _update_results_step(self):
        ipw.dlink(
            (self.submit_model, "process_node"),
            (self.results_model, "process_uuid"),
            lambda node: node.uuid if node is not None else None,
        )

    def _lock_app(self):
        for model in (
            self.structure_model,
            self.configure_model,
            self.submit_model,
        ):
            model.unobserve_all("confirmed")

    def _update_from_process(self, pk):
        if pk is None:
            self._wizard_app_widget.reset()
            self._wizard_app_widget.selected_index = 0
        else:
            self._show_process_loading_message()
            process_node = load_node(pk)
            self.structure_model.input_structure = process_node.inputs.structure
            self.structure_model.confirm()
            parameters = process_node.base.extras.get("ui_parameters", {})
            if parameters and isinstance(parameters, str):
                parameters = deserialize_unsafe(parameters)
            self.configure_model.set_model_state(parameters)
            self.configure_model.confirm()
            self.submit_model.process_node = process_node
            self.submit_model.set_model_state(parameters)
            self.submit_model.confirm()
            self._wizard_app_widget.selected_index = 3
            self._hide_process_loading_message()

    def _show_process_loading_message(self):
        self._process_loading_message.layout.display = "flex"

    def _hide_process_loading_message(self):
        self._process_loading_message.layout.display = "none"
