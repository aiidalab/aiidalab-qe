"""Class to .

Authors:

    AiiDAlab Team
"""

from __future__ import annotations

import os
import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.code.model import CodeModel
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.mixins import Confirmable, HasBlockers, HasModels, HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import (
    LoadingWidget,
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class PanelModel(Model):
    """Base class for all panel models.

    Attributes
    ----------
    `title` : `str`
        The title to be shown in the GUI.
    `identifier` : `str`
        Which plugin this panel belong to.
    """

    title = ""
    identifier = ""


PM = t.TypeVar("PM", bound=PanelModel)


class Panel(ipw.VBox, t.Generic[PM]):
    """Base class for all panels."""

    rendered = False
    loading_message = "Loading {identifier} panel"

    def __init__(self, model: PM, **kwargs):
        loading_message = self.loading_message.format(identifier=model.identifier)
        loading_message = loading_message.replace("_", " ")
        self.loading_message = LoadingWidget(loading_message)
        super().__init__(children=[self.loading_message], **kwargs)
        self._model = model

    def render(self):
        raise NotImplementedError()


class PluginOutline(ipw.HBox):
    title = "Outline"
    description = ""

    def __init__(self, **kwargs):
        self.include = ipw.Checkbox(
            description=self.title,
            indent=False,
            style={"description_width": "initial"},
        )

        super().__init__(
            children=[
                self.include,
                ipw.HTML(f"""
                    <div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
                        {self.description}
                    </div>
                """),
            ],
            **kwargs,
        )


class SettingsModel(PanelModel, HasBlockers):
    """Base model for settings models."""

    dependencies: list[str] = []

    include = tl.Bool(False)
    loaded_from_process = tl.Bool(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._defaults = {}

    def update(self):
        """Updates the model."""
        pass

    def get_model_state(self) -> dict:
        """Retrieves the current state of the model as a dictionary."""
        raise NotImplementedError()

    def set_model_state(self, parameters: dict):
        """Distributes the parameters of a loaded calculation to the model."""
        raise NotImplementedError()

    def reset(self):
        """Resets the model to present defaults."""
        pass


SM = t.TypeVar("SM", bound=SettingsModel)


class SettingsPanel(Panel[SM]):
    """Base model for settings panels."""

    updated = False

    def __init__(self, model: SM, **kwargs):
        super().__init__(model=model, **kwargs)
        self.links = []


class ConfigurationSettingsModel(SettingsModel, Confirmable):
    """Base model for configuration settings models."""

    def update(self, specific=""):
        """Updates the model.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of update.
        """
        pass


CSM = t.TypeVar("CSM", bound=ConfigurationSettingsModel)


class ConfigurationSettingsPanel(SettingsPanel[CSM]):
    """Base class for configuration settings panels."""

    def refresh(self, specific=""):
        """Refreshes the settings panel.

        Unlinks the panel's widgets. If the panel's model is included in the
        calculation, also updates the model's defaults. Resets the model to
        these defaults if there is no input structure.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of refresh.
        """
        self.updated = False
        self._unsubscribe()
        if self._model.include:
            self.update(specific)
        if "PYTEST_CURRENT_TEST" in os.environ:
            # Skip resetting to avoid having to inject a structure when testing
            return
        if hasattr(self._model, "input_structure") and not self._model.input_structure:
            self._reset()

    def update(self, specific=""):
        """Updates the model if not yet updated.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of update.
        """
        if self.updated:
            return
        if not self._model.loaded_from_process:
            self._model.update(specific)
        self.updated = True

    def _unsubscribe(self):
        """Unlinks any linked widgets."""
        for link in self.links:
            link.unlink()
        self.links.clear()

    def _reset(self):
        """Resets the model to present defaults."""
        self.updated = False
        self._model.reset()


class ResourceSettingsModel(SettingsModel, HasModels[CodeModel]):
    """Base model for resource setting models."""

    global_codes = tl.Dict(
        key_trait=tl.Unicode(),
        value_trait=tl.Dict(),
    )

    warning_messages = tl.Unicode("")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Used by the code-setup thread to fetch code options
        self.DEFAULT_USER_EMAIL = orm.User.collection.get_default().email

    def add_model(self, identifier, model):
        super().add_model(identifier, model)
        model.update(self.DEFAULT_USER_EMAIL)

    def refresh_codes(self):
        for _, code_model in self.get_models():
            code_model.update(self.DEFAULT_USER_EMAIL, refresh=True)

    def get_model_state(self):
        return {
            "codes": {
                identifier: code_model.get_model_state()
                for identifier, code_model in self.get_models()
                if code_model.is_ready
            },
        }

    def set_model_state(self, parameters: dict):
        self.set_selected_codes(parameters.get("codes", {}))

    def get_selected_codes(self) -> dict[str, dict]:
        return {
            identifier: code_model.get_model_state()
            for identifier, code_model in self.get_models()
            if code_model.is_ready
        }

    def set_selected_codes(self, code_data=DEFAULT["codes"]):
        for identifier, code_model in self.get_models():
            if identifier in code_data:
                code_model.set_model_state(code_data[identifier])

    def _check_blockers(self):
        return []


RSM = t.TypeVar("RSM", bound=ResourceSettingsModel)


class ResourceSettingsPanel(SettingsPanel[RSM]):
    """Base class for resource setting panels."""

    def __init__(self, model, **kwargs):
        super().__init__(model, **kwargs)
        self.code_widgets = {}

    def _on_code_resource_change(self, _):
        pass

    def _on_code_options_change(self, change: dict):
        widget: ipw.Dropdown = change["owner"]
        widget.disabled = not widget.options

    def _toggle_code(self, code_model: CodeModel):
        if not self.rendered:
            return
        if not code_model.is_rendered:
            loading_message = LoadingWidget(f"Loading {code_model.name} code")
            self.code_widgets_container.children += (loading_message,)
        if code_model.name not in self.code_widgets:
            code_widget = code_model.code_widget_class(
                description=code_model.description,
                default_calc_job_plugin=code_model.default_calc_job_plugin,
            )
            self.code_widgets[code_model.name] = code_widget
        else:
            code_widget = self.code_widgets[code_model.name]
        if not code_model.is_rendered:
            code_widget.observe(
                code_widget.update_resources,
                "value",
            )
            self._render_code_widget(code_model, code_widget)

    def _render_code_widget(
        self,
        code_model: CodeModel,
        code_widget: QEAppComputationalResourcesWidget,
    ):
        ipw.dlink(
            (code_model, "options"),
            (code_widget.code_selection.code_select_dropdown, "options"),
        )
        ipw.link(
            (code_model, "selected"),
            (code_widget, "value"),
        )
        ipw.link(
            (code_model, "num_cpus"),
            (code_widget.num_cpus, "value"),
        )
        ipw.link(
            (code_model, "num_nodes"),
            (code_widget.num_nodes, "value"),
        )
        ipw.link(
            (code_model, "ntasks_per_node"),
            (code_widget.resource_detail.ntasks_per_node, "value"),
        )
        ipw.link(
            (code_model, "cpus_per_task"),
            (code_widget.resource_detail.cpus_per_task, "value"),
        )
        ipw.link(
            (code_model, "max_wallclock_seconds"),
            (code_widget.resource_detail.max_wallclock_seconds, "value"),
        )
        if isinstance(code_widget, PwCodeResourceSetupWidget):
            ipw.link(
                (code_model, "parallelization_override"),
                (code_widget.parallelization.override, "value"),
            )
            ipw.link(
                (code_model, "npool"),
                (code_widget.parallelization.npool, "value"),
            )
            code_model.observe(
                self._on_code_resource_change,
                [
                    "parallelization_override",
                    "npool",
                ],
            )
        code_model.observe(
            self._on_code_resource_change,
            [
                "selected",
                "num_cpus",
                "num_nodes",
                "ntasks_per_node",
                "cpus_per_task",
                "max_wallclock_seconds",
            ],
        )
        code_widget.code_selection.code_select_dropdown.observe(
            self._on_code_options_change,
            "options",
        )
        code_widgets = self.code_widgets_container.children[:-1]  # type: ignore
        self.code_widgets_container.children = [*code_widgets, code_widget]
        code_model.is_rendered = True


class PluginResourceSettingsModel(ResourceSettingsModel):
    """Base model for plugin resource setting models."""

    dependencies = [
        "global.global_codes",
    ]

    override = tl.Bool(False)

    def add_model(self, identifier, model: CodeModel):
        super().add_model(identifier, model)
        model.activate()

    def update(self):
        """Updates the code models from the global resources.

        Skips synchronization with global resources if the user has chosen to override
        the resources for the plugin codes.
        """
        if self.override:
            return
        for _, code_model in self.get_models():
            model_key = code_model.default_calc_job_plugin.replace(".", "__")
            if model_key in self.global_codes:
                code_resources: dict = self.global_codes[model_key]  # type: ignore
                code_model.set_model_state(code_resources)

    def get_model_state(self):
        return {
            "override": self.override,
            **super().get_model_state(),
        }

    def set_model_state(self, parameters: dict):
        self.override = parameters.get("override", False)
        super().set_model_state(parameters)

    def _link_model(self, model: CodeModel):
        tl.link(
            (self, "override"),
            (model, "override"),
        )


PRSM = t.TypeVar("PRSM", bound=PluginResourceSettingsModel)


class PluginResourceSettingsPanel(ResourceSettingsPanel[PRSM]):
    """Base class for plugin resource setting panels."""

    def __init__(self, model, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_global_codes_change,
            "global_codes",
        )
        self._model.observe(
            self._on_override_change,
            "override",
        )

    def render(self):
        if self.rendered:
            return

        self.override_help = ipw.HTML(
            "Click to override the resource settings for this plugin."
        )
        self.override = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="3%"),
        )
        ipw.link(
            (self._model, "override"),
            (self.override, "value"),
        )
        self.code_widgets_container = ipw.VBox()

        self.children = [
            ipw.HBox(
                children=[
                    self.override,
                    self.override_help,
                ]
            ),
            self.code_widgets_container,
        ]

        self.rendered = True

        # Render any active codes
        for _, code_model in self._model.get_models():
            if code_model.is_active:
                self._toggle_code(code_model)

        return self.code_widgets_container

    def _on_global_codes_change(self, _):
        self._model.update()

    def _on_override_change(self, _):
        self._model.update()

    def _render_code_widget(
        self,
        code_model: CodeModel,
        code_widget: QEAppComputationalResourcesWidget,
    ):
        super()._render_code_widget(code_model, code_widget)
        self._link_override_to_widget_disable(code_model, code_widget)

    def _link_override_to_widget_disable(self, code_model, code_widget):
        """Links the override attribute of the code model to the disable attribute
        of subwidgets of the code widget."""
        ipw.dlink(
            (code_model, "override"),
            (code_widget.code_selection.code_select_dropdown, "disabled"),
            lambda override: not override,
        )
        ipw.dlink(
            (code_model, "override"),
            (code_widget.num_cpus, "disabled"),
            lambda override: not override,
        )
        ipw.dlink(
            (code_model, "override"),
            (code_widget.num_nodes, "disabled"),
            lambda override: not override,
        )
        ipw.dlink(
            (code_model, "override"),
            (code_widget.btn_setup_resource_detail, "disabled"),
            lambda override: not override,
        )
        if isinstance(code_widget, PwCodeResourceSetupWidget):
            ipw.dlink(
                (code_model, "override"),
                (code_widget.parallelization.override, "disabled"),
                lambda override: not override,
            )
            ipw.dlink(
                (code_model, "override"),
                (code_widget.parallelization.npool, "disabled"),
                lambda override: not override,
            )


class ResultsModel(PanelModel, HasProcess):
    process_status_notification = tl.Unicode("")

    _this_process_label = ""
    _this_process_uuid = None

    auto_render = False
    _completed_process = False

    CSS_MAP = {
        "finished": "success",
        "failed": "danger",
        "excepted": "danger",
        "killed": "danger",
        "queued": "warning",
        "running": "info",
        "created": "info",
    }

    @property
    def include(self):
        return self.identifier in self.properties

    @property
    def has_results(self):
        node = self.fetch_child_process_node()
        return node and node.is_finished_ok

    def update(self):
        self.auto_render = self.has_results

    def update_process_status_notification(self):
        if self._completed_process:
            self.process_status_notification = ""
            return
        status = self._get_child_process_status()
        self.process_status_notification = status
        if "success" in status:
            self._completed_process = True

    def fetch_child_process_node(self, which="this") -> orm.ProcessNode | None:
        if not self.process_uuid:
            return
        which = which.lower()
        uuid = getattr(self, f"_{which}_process_uuid")
        label = getattr(self, f"_{which}_process_label")
        if not uuid:
            root = self.fetch_process_node()
            child = next((c for c in root.called if c.process_label == label), None)
            uuid = child.uuid if child else None
        return orm.load_node(uuid) if uuid else None  # type: ignore

    def _get_child_process_status(self, which="this"):
        state, exit_message = self._get_child_state_and_exit_message(which)
        if state == "waiting":
            state = "running"
        status = state.upper()
        if exit_message:
            status = f"{status} ({exit_message})"
        label = "Status" if which == "this" else f"{which.capitalize()} status"
        alert_class = f"alert-{self.CSS_MAP.get(state, 'info')}"
        return f"""
            <div class="alert {alert_class}" style="padding: 5px 10px;">
                <b>{label}:</b> {status}
            </div>
        """

    def _get_child_state_and_exit_message(self, which="this"):
        if not (
            (node := self.fetch_child_process_node(which))
            and hasattr(node, "process_state")
            and node.process_state
        ):
            return "queued", None
        if node.is_failed:
            return "failed", node.exit_message
        return node.process_state.value, None

    def _get_child_outputs(self, which="this"):
        if not (node := self.fetch_child_process_node(which)):
            outputs = super().outputs
            child = which if which != "this" else self.identifier
            return getattr(outputs, child) if child in outputs else AttributeDict({})
        return AttributeDict({key: getattr(node.outputs, key) for key in node.outputs})


RM = t.TypeVar("RM", bound=ResultsModel)


class ResultsPanel(Panel[RM]):
    """Base class for all the result panels.

    The base class has a method to load the result of the calculation.
    And a show method to display it in the panel.
    It has a update method to update the result in the panel.
    """

    loading_message = "Loading {identifier} results"

    def __init__(self, model: RM, **kwargs):
        super().__init__(model=model, **kwargs)
        self._model.observe(
            self._on_process_change,
            "process_uuid",
        )
        self._model.observe(
            self._on_monitor_counter_change,
            "monitor_counter",
        )

    def render(self):
        if self.rendered:
            if self._model.identifier == "structure":
                self._render()
            return

        if not self._model.has_process:
            return

        self.guide = InAppGuide(
            identifier=f"{self._model.identifier}-results",
            classes=["results-panel-guide"],
        )

        self.results_container = ipw.VBox()

        if self._model.auto_render:
            self.children = [
                self.guide,
                self.results_container,
            ]
            self._load_results()
        else:
            children = [self.guide]
            if (
                self._model.identifier != "structure"
                or "relax" in self._model.properties
            ):
                children.append(self._get_controls_section())
            children.append(self.results_container)
            self.children = children
            if self._model.identifier == "structure":
                self._load_results()

        self.rendered = True

    def _on_process_change(self, _):
        self._model.update()

    def _on_monitor_counter_change(self, _):
        self._model.update_process_status_notification()

    def _on_load_results_click(self, _):
        self.load_controls.children = []
        self._load_results()

    def _load_results(self):
        self.results_container.children = [self.loading_message]
        self._render()
        self._post_render()

    def _get_controls_section(self) -> ipw.VBox:
        self.process_status_notification = ipw.HTML()
        ipw.dlink(
            (self._model, "process_status_notification"),
            (self.process_status_notification, "value"),
        )

        self.load_results_button = ipw.Button(
            description="Load results",
            button_style="warning",
            tooltip="Load the results",
            icon="refresh",
        )
        ipw.dlink(
            (self._model, "monitor_counter"),
            (self.load_results_button, "disabled"),
            lambda _: not self._model.has_results,
        )
        self.load_results_button.on_click(self._on_load_results_click)

        self.load_controls = ipw.HBox(
            children=[]
            if self._model.auto_render or self._model.identifier == "structure"
            else [
                self.load_results_button,
                ipw.HTML("""
                    <div style="margin-left: 10px">
                        <b>Note:</b> Load time may vary depending on the size of the
                        calculation
                    </div>
                """),
            ]
        )

        return ipw.VBox(
            children=[
                self.process_status_notification,
                self.load_controls,
            ]
        )

    def _render(self):
        raise NotImplementedError()

    def _post_render(self):
        pass
