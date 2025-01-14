"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.code import PluginCodes, PwCodeModel
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import (
    PluginResourceSettingsModel,
    PluginResourceSettingsPanel,
    ResourceSettingsPanel,
)
from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.common.setup_pseudos import PseudosInstallWidget
from aiidalab_qe.common.widgets import LinkButton, QeDependentWizardStep

from .global_settings import GlobalResourceSettingsModel, GlobalResourceSettingsPanel
from .model import SubmissionStepModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmitQeAppWorkChainStep(QeDependentWizardStep[SubmissionStepModel]):
    missing_information_warning = "Missing input structure and/or configuration parameters. Please set them first."

    def __init__(self, model: SubmissionStepModel, qe_auto_setup=True, **kwargs):
        super().__init__(model=model, **kwargs)
        self._model.observe(
            self._on_submission,
            "confirmed",
        )
        self._model.observe(
            self._on_input_parameters_change,
            "input_parameters",
        )
        self._model.observe(
            self._on_submission_blockers_change,
            [
                "internal_submission_blockers",
                "external_submission_blockers",
            ],
        )
        self._model.observe(
            self._on_installation_change,
            ["installing_sssp", "sssp_installed"],
        )
        self._model.observe(
            self._on_sssp_installed,
            "sssp_installed",
        )
        self._model.observe(
            self._on_installation_change,
            ["installing_qe", "qe_installed"],
        )
        self._model.observe(
            self._on_qe_installed,
            "qe_installed",
        )

        global_resources_model = GlobalResourceSettingsModel()
        self.global_resources = GlobalResourceSettingsPanel(
            model=global_resources_model
        )
        self._model.add_model("global", global_resources_model)
        ipw.dlink(
            (self._model, "plugin_overrides"),
            (global_resources_model, "plugin_overrides"),
        )
        global_resources_model.observe(
            self._on_plugin_submission_blockers_change,
            ["submission_blockers"],
        )
        global_resources_model.observe(
            self._on_plugin_submission_warning_messages_change,
            ["submission_warning_messages"],
        )

        self.settings = {
            "global": self.global_resources,
        }
        self._fetch_plugin_resource_settings()

        self._install_sssp(qe_auto_setup)
        self._set_up_qe(qe_auto_setup)

    def _render(self):
        self.process_label = ipw.Text(
            description="Label:",
            layout=ipw.Layout(width="auto", indent="0px"),
        )
        ipw.link(
            (self._model, "process_label"),
            (self.process_label, "value"),
        )
        self.process_description = ipw.Textarea(
            description="Description",
            layout=ipw.Layout(width="auto", indent="0px"),
        )
        ipw.link(
            (self._model, "process_description"),
            (self.process_description, "value"),
        )

        self.submit_button = ipw.Button(
            description="Submit",
            tooltip="Submit the calculation with the selected parameters.",
            icon="play",
            button_style="success",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )
        ipw.dlink(
            (self, "state"),
            (self.submit_button, "disabled"),
            lambda state: state != self.State.CONFIGURED,
        )
        self.submit_button.on_click(self.submit)

        self.submission_blocker_messages = ipw.HTML()
        ipw.dlink(
            (self._model, "submission_blocker_messages"),
            (self.submission_blocker_messages, "value"),
        )

        self.submission_warning_messages = ipw.HTML()
        ipw.dlink(
            (self._model, "submission_warning_messages"),
            (self.submission_warning_messages, "value"),
        )

        self.setup_new_codes_button = LinkButton(
            description="Setup resources",
            link="../home/code_setup.ipynb",
            icon="database",
        )

        self.refresh_resources_button = ipw.Button(
            description="Refresh resources",
            icon="refresh",
            button_style="primary",
            layout=ipw.Layout(width="fit-content", margin="2px 2px 12px"),
        )
        self.refresh_resources_button.on_click(self._refresh_resources)

        self.tabs = ipw.Tab(
            layout=ipw.Layout(min_height="250px"),
            selected_index=None,
        )
        self.tabs.observe(
            self._on_tab_change,
            "selected_index",
        )

        self.children = [
            InAppGuide(identifier="submission-step"),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
                    Select the codes to use for running the calculations. The codes on
                    the local machine (localhost) are installed by default, but you can
                    configure new ones on potentially more powerful machines by clicking
                    on <i class="fa fa-database"></i> <b>Setup resources</b> (also at
                    the top of the app). Make sure to click the <b>Refresh resources</b>
                    button below after making changes to AiiDA resources to update the
                    app resources.
                </div>
            """),
            ipw.HBox(
                children=[
                    self.setup_new_codes_button,
                    self.refresh_resources_button,
                ],
                layout=ipw.Layout(grid_gap="5px"),
            ),
            self.tabs,
            self.sssp_installation,
            self.qe_setup,
            self.submission_blocker_messages,
            self.submission_warning_messages,
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
                    Label your job and provide a brief description. These details
                    help identify the job later and make the search process easier.
                    While optional, adding a description is recommended for better
                    clarity.
                </div>
            """),
            self.process_label,
            self.process_description,
            self.submit_button,
        ]

    def _post_render(self):
        self._update_tabs()

    def submit(self, _=None):
        self._model.confirm()

    def reset(self):
        self._model.reset()

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update_state()

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: ResourceSettingsPanel = self.tabs.children[tab_index]  # type: ignore
        tab.render()

    def _on_input_parameters_change(self, _):
        self._model.update_process_label()
        self._model.update_plugin_inclusion()
        self._model.update_plugin_overrides()
        self._model.update_submission_blockers()
        self._update_tabs()

    def _on_plugin_overrides_change(self, _):
        self._model.update_plugin_overrides()

    def _on_plugin_submission_blockers_change(self, _):
        self._model.update_submission_blockers()

    def _on_plugin_submission_warning_messages_change(self, _):
        self._model.update_submission_warnings()

    def _on_submission_blockers_change(self, _):
        self._model.update_submission_blocker_message()
        self._update_state()

    def _on_installation_change(self, _):
        self._model.update_submission_blockers()

    def _on_qe_installed(self, _):
        self._toggle_qe_installation_widget()
        if self._model.qe_installed:
            self._model.update()

    def _on_sssp_installed(self, _):
        self._toggle_sssp_installation_widget()

    def _on_submission(self, _):
        self._update_state()

    def _install_sssp(self, qe_auto_setup):
        self.sssp_installation = PseudosInstallWidget(auto_start=False)
        ipw.dlink(
            (self.sssp_installation, "busy"),
            (self._model, "installing_sssp"),
        )
        ipw.dlink(
            (self.sssp_installation, "installed"),
            (self._model, "installing_sssp"),
            lambda installed: not installed,
        )
        ipw.dlink(
            (self.sssp_installation, "installed"),
            (self._model, "sssp_installed"),
        )
        if qe_auto_setup:
            self.sssp_installation.refresh()

    def _set_up_qe(self, qe_auto_setup):
        self.qe_setup = QESetupWidget(auto_start=False)
        ipw.dlink(
            (self.qe_setup, "busy"),
            (self._model, "installing_qe"),
        )
        ipw.dlink(
            (self.qe_setup, "installed"),
            (self._model, "installing_qe"),
            lambda installed: not installed,
        )
        ipw.dlink(
            (self.qe_setup, "installed"),
            (self._model, "qe_installed"),
        )
        if qe_auto_setup:
            self.qe_setup.refresh()

    def _toggle_sssp_installation_widget(self):
        sssp_installation_display = "none" if self._model.sssp_installed else "block"
        self.sssp_installation.layout.display = sssp_installation_display

    def _toggle_qe_installation_widget(self):
        qe_installation_display = "none" if self._model.qe_installed else "block"
        self.qe_setup.layout.display = qe_installation_display

    def _refresh_resources(self, _=None):
        for _, model in self._model.get_models():
            model.refresh_codes()

    def _update_tabs(self):
        children = []
        titles = []
        for identifier, model in self._model.get_models():
            if model.include:
                settings = self.settings[identifier]
                titles.append(model.title)
                children.append(settings)
        if self.rendered:
            self.tabs.selected_index = None
            self.tabs.children = children
            for i, title in enumerate(titles):
                self.tabs.set_title(i, title)
            self.tabs.selected_index = 0

    def _update_state(self, _=None):
        if self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
        elif self.previous_step_state is not self.State.SUCCESS:
            self.state = self.State.INIT
        elif self._model.confirmed:
            self.state = self.State.SUCCESS
        elif self._model.is_blocked:
            self.state = self.State.READY
        else:
            self.state = self.state.CONFIGURED

    def _fetch_plugin_resource_settings(self):
        entries = get_entry_items("aiidalab_qe.properties", "resources")
        codes: PluginCodes = {
            "dft": {
                "pw": PwCodeModel(),
            },
        }
        for identifier, resources in entries.items():
            for key in ("panel", "model"):
                if key not in resources:
                    raise ValueError(f"Entry {identifier} is missing the '{key}' key")

            model: PluginResourceSettingsModel = resources["model"]()
            model.observe(
                self._on_plugin_overrides_change,
                "override",
            )
            model.observe(
                self._on_plugin_submission_blockers_change,
                ["submission_blockers"],
            )
            model.observe(
                self._on_plugin_submission_warning_messages_change,
                ["submission_warning_messages"],
            )
            self._model.add_model(identifier, model)

            panel: PluginResourceSettingsPanel = resources["panel"](model=model)
            self.settings[identifier] = panel

            codes[identifier] = dict(model.get_models())

        self.global_resources.set_up_codes(codes)
