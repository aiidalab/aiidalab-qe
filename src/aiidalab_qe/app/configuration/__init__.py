"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.panel import SettingsModel, SettingsPanel
from aiidalab_widgets_base import WizardAppWidgetStep

from .advanced import AdvancedModel, AdvancedSettings
from .basic import WorkChainModel, WorkChainSettings
from .model import ConfigurationModel

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class ConfigureQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    def __init__(self, model: ConfigurationModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget("Loading workflow configuration panel")],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_confirmation_change,
            "confirmed",
        )
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

        self.structure_set_message = ipw.HTML()
        ipw.dlink(
            (self._model, "input_structure"),
            (self.structure_set_message, "value"),
            lambda structure: ""
            if structure
            else """
                <div class="alert alert-info">
                    <b>Please set the input structure first.</b>
                </div>
            """,
        )

        workchain_model = WorkChainModel()
        self.workchain_settings = WorkChainSettings(model=workchain_model)
        self._model.add_model("workchain", workchain_model)

        advanced_model = AdvancedModel()
        self.advanced_settings = AdvancedSettings(model=advanced_model)
        self._model.add_model("advanced", advanced_model)

        self.settings = {
            "workchain": self.workchain_settings,
            "advanced": self.advanced_settings,
        }

        self._fetch_plugin_settings()

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type_help = ipw.HTML()
        ipw.dlink(
            (self._model, "relax_type_help"),
            (self.relax_type_help, "value"),
        )
        self.relax_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "none"),
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ],
        )
        ipw.dlink(
            (self._model, "relax_type_options"),
            (self.relax_type, "options"),
        )
        ipw.link(
            (self._model, "relax_type"),
            (self.relax_type, "value"),
        )

        self.tab = ipw.Tab(
            layout=ipw.Layout(min_height="250px"),
            selected_index=None,
        )
        self.tab.observe(
            self._on_tab_change,
            "selected_index",
        )
        self._update_tabs()

        self.sub_steps = ipw.Accordion(
            children=[
                ipw.VBox(
                    children=[
                        *self.property_children,
                    ]
                ),
                self.tab,
            ],
            layout=ipw.Layout(margin="10px 2px"),
            selected_index=None,
        )
        self.sub_steps.set_title(0, "Step 2.1: Select which properties to calculate")
        self.sub_steps.set_title(1, "Step 2.2: Customize calculation parameters")

        self.confirm_button = ipw.Button(
            description="Confirm",
            tooltip="Confirm the currently selected settings and go to the next step.",
            button_style="success",
            icon="check-circle",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        ipw.dlink(
            (self, "state"),
            (self.confirm_button, "disabled"),
            lambda state: state != self.State.CONFIGURED,
        )
        self.confirm_button.on_click(self.confirm)

        self.children = [
            self.structure_set_message,
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Structure relaxation</h4>
                </div>
            """),
            self.relax_type_help,
            self.relax_type,
            self.sub_steps,
            self.confirm_button,
        ]

        self.rendered = True

        self._model.update()

    def is_saved(self):
        return self._model.confirmed

    def confirm(self, _=None):
        self._model.confirm()

    def reset(self):
        self._model.reset()
        if self.rendered:
            self.sub_steps.selected_index = None
            self.tab.selected_index = 0

    @tl.observe("previous_step_state")
    def _on_previous_step_state_change(self, _):
        self._update_state()

    def _on_tab_change(self, change):
        if (tab_index := change["new"]) is None:
            return
        tab: SettingsPanel = self.tab.children[tab_index]  # type: ignore
        tab.render()
        tab.update()

    def _on_input_structure_change(self, _):
        self._model.update()
        self.reset()

    def _on_confirmation_change(self, _):
        self._update_state()

    def _update_tabs(self):
        children = []
        titles = []
        for identifier, model in self._model.get_models():
            if model.include:
                settings = self.settings[identifier]
                titles.append(settings.title)
                children.append(settings)
        if hasattr(self, "tab"):
            self.tab.selected_index = None
            self.tab.children = children
            for i, title in enumerate(titles):
                self.tab.set_title(i, title)
            self.tab.selected_index = 0

    def _update_state(self, _=None):
        if self._model.confirmed:
            self.state = self.State.SUCCESS
        elif self.previous_step_state is self.State.SUCCESS:
            self.state = self.State.CONFIGURED
        elif self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
        else:
            self.state = self.State.INIT

    def _fetch_plugin_settings(self):
        self.property_children = []

        outlines = get_entry_items("aiidalab_qe.properties", "outline")
        models = get_entry_items("aiidalab_qe.properties", "model")
        settings = get_entry_items("aiidalab_qe.properties", "setting")
        for identifier in settings:
            model: SettingsModel = models[identifier]()
            self._model.add_model(identifier, model)

            outline = outlines[identifier]()
            info = ipw.HTML()
            ipw.link(
                (model, "include"),
                (outline.include, "value"),
            )

            if identifier == "bands":
                ipw.dlink(
                    (self._model, "input_structure"),
                    (outline.include, "disabled"),
                    lambda _: not self._model.has_pbc,
                )

            def toggle_plugin(change, identifier=identifier, model=model, info=info):
                if change["new"]:
                    info.value = (
                        f"Customize {identifier} settings in <b>Step 2.2</b> if needed"
                    )
                else:
                    info.value = ""
                model.update()
                self._update_tabs()

            model.observe(
                toggle_plugin,
                "include",
            )

            self.property_children.append(
                ipw.HBox(
                    children=[
                        outline,
                        info,
                    ]
                )
            )

            self.settings[identifier] = settings[identifier](
                parent=self,
                identifier=identifier,
                model=model,
            )
