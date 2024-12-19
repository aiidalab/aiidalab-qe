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
from aiidalab_qe.common.mixins import Confirmable, HasModels, HasProcess
from aiidalab_qe.common.mvc import Model
from aiidalab_qe.common.widgets import (
    LoadingWidget,
    PwCodeResourceSetupWidget,
    QEAppComputationalResourcesWidget,
)

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class Panel(ipw.VBox):
    """Base class for all the panels.

    The base class has a method to return the value of all the widgets in
    the panel as a dictionary. The dictionary is used to construct the
    input file for the calculation. The class also has a method to load a dictionary to set the value of the widgets in the panel.

    title: the title to be shown in the GUI
    identifier: which plugin this panel belong to.

    """

    title = "Panel"

    # TODO remove `identifier` (and `parent`) from signature
    # TODO add `model` parameter
    # TODO add `identifier` property to route to model.identifier
    def __init__(self, parent=None, identifier=None, **kwargs):
        """Initialize the panel.

        :param kwargs: keyword arguments to pass to the ipw.VBox constructor.
        """
        self.parent = parent
        self.identifier = identifier or getattr(self, "identifier", "plugin")
        super().__init__(
            children=kwargs.pop("children", []),
            **kwargs,
        )


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


class SettingsModel(Model):
    title = "Model"
    identifier = ""
    dependencies: list[str] = []

    include = tl.Bool(False)
    loaded_from_process = tl.Bool(False)

    _defaults = {}

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


class SettingsPanel(Panel, t.Generic[SM]):
    title = "Settings"

    def __init__(self, model: SM, **kwargs):
        self.loading_message = LoadingWidget(f"Loading {model.identifier} settings")

        super().__init__(
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model

        self.rendered = False
        self.updated = False

        self.links = []

    def render(self):
        raise NotImplementedError()


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


class ConfigurationSettingsPanel(SettingsPanel[CSM], t.Generic[CSM]):
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

    dependencies = []

    global_codes = tl.Dict(
        key_trait=tl.Unicode(),
        value_trait=tl.Dict(),
    )

    submission_blockers = tl.List(tl.Unicode())
    submission_warning_messages = tl.Unicode("")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Used by the code-setup thread to fetch code options
        self.DEFAULT_USER_EMAIL = orm.User.collection.get_default().email

    def add_model(self, identifier, model):
        super().add_model(identifier, model)
        model.update(self.DEFAULT_USER_EMAIL)

    def update_submission_blockers(self):
        self.submission_blockers = list(self._check_submission_blockers())

    def get_model_state(self):
        return {
            "codes": {
                identifier: code_model.get_model_state()
                for identifier, code_model in self.get_models()
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

    def _check_submission_blockers(self):
        return []


RSM = t.TypeVar("RSM", bound=ResourceSettingsModel)


class ResourceSettingsPanel(SettingsPanel[RSM], t.Generic[RSM]):
    """Base class for resource setting panels."""

    def __init__(self, model, **kwargs):
        super().__init__(model, **kwargs)

        self.code_widgets = {}

    def _on_code_resource_change(self, _):
        pass

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
            self._render_code_widget(code_model, code_widget)
            code_widget.observe(
                code_widget.update_resources,
                "value",
            )

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
                "options",
                "selected",
                "num_cpus",
                "num_nodes",
                "ntasks_per_node",
                "cpus_per_task",
                "max_wallclock_seconds",
            ],
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

    def update(self):
        """Updates the code models from the global resources.

        Skips synchronization with global resources if the user has chosen to override
        the resources for the plugin codes.
        """
        if self.override:
            return
        for _, code_model in self.get_models():
            default_calc_job_plugin = code_model.default_calc_job_plugin
            if default_calc_job_plugin in self.global_codes:
                code_resources: dict = self.global_codes[default_calc_job_plugin]  # type: ignore
                options = code_resources.get("options", [])
                if options != code_model.options:
                    code_model.update(self.DEFAULT_USER_EMAIL, refresh=True)
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


class PluginResourceSettingsPanel(ResourceSettingsPanel[PRSM], t.Generic[PRSM]):
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
            (code_widget.code_selection.btn_setup_new_code, "disabled"),
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


class ResultsModel(Model, HasProcess):
    title = "Model"
    identifier = "model"

    process_status_notification = tl.Unicode("")

    _this_process_label = ""
    _this_process_uuid = None

    CSS_MAP = {
        "finished": "success",
        "failed": "danger",
        "excepted": "danger",
        "killed": "danger",
        "queued": "warning",
        "waiting": "info",
        "running": "info",
        "created": "info",
    }

    @property
    def has_results(self):
        node = self._fetch_child_process_node()
        return node and node.is_finished_ok

    def update_process_status_notification(self):
        self.process_status_notification = self._get_child_process_status()

    def _get_child_process_status(self, child="this"):
        state, exit_message = self._get_child_state_and_exit_message(child)
        status = state.upper()
        if exit_message:
            status = f"{status} ({exit_message})"
        label = "Status" if child == "this" else f"{child.capitalize()} status"
        alert_class = f"alert-{self.CSS_MAP.get(state, 'info')}"
        return f"""
            <div class="alert {alert_class}" style="padding: 5px 10px;">
                <b>{label}:</b> {status}
            </div>
        """

    def _get_child_state_and_exit_message(self, child="this"):
        if not (
            (node := self._fetch_child_process_node(child))
            and hasattr(node, "process_state")
            and node.process_state
        ):
            return "queued", None
        if node.is_failed:
            return "failed", node.exit_message
        return node.process_state.value, None

    def _get_child_outputs(self, child="this"):
        if not (node := self._fetch_child_process_node(child)):
            outputs = super().outputs
            child = child if child != "this" else self.identifier
            return getattr(outputs, child) if child in outputs else AttributeDict({})
        return AttributeDict({key: getattr(node.outputs, key) for key in node.outputs})

    def _fetch_child_process_node(self, child="this") -> orm.ProcessNode | None:
        if not self.process_uuid:
            return
        child = child.lower()
        uuid = getattr(self, f"_{child}_process_uuid")
        label = getattr(self, f"_{child}_process_label")
        if not uuid:
            uuid = (
                orm.QueryBuilder()
                .append(
                    orm.WorkChainNode,
                    filters={"uuid": self.process_uuid},
                    tag="root_process",
                )
                .append(
                    orm.WorkChainNode,
                    filters={"attributes.process_label": label},
                    project="uuid",
                    with_incoming="root_process",
                )
                .first(flat=True)
            )
        return orm.load_node(uuid) if uuid else None  # type: ignore


RM = t.TypeVar("RM", bound=ResultsModel)


class ResultsPanel(Panel, t.Generic[RM]):
    """Base class for all the result panels.

    The base class has a method to load the result of the calculation.
    And a show method to display it in the panel.
    It has a update method to update the result in the panel.
    """

    title = "Results"
    identifier = "results"

    # To specify which plugins (outputs) are needed
    # for this result panel.
    workchain_labels = []

    def __init__(self, model: RM, **kwargs):
        self.loading_message = LoadingWidget(f"Loading {self.title.lower()} results")

        super().__init__(**kwargs)

        self._model = model
        self._model.observe(
            self._on_process_change,
            "process_uuid",
        )
        self._model.observe(
            self._on_monitor_counter_change,
            "monitor_counter",
        )

        self.rendered = False
        self.has_controls = False

        self.links = []

    def render(self):
        if self.rendered:
            if self.identifier == "structure":
                self._render()
            return
        if self.has_controls or not self._model.has_process:
            return
        if not self._model.has_results:
            self._render_controls()
        else:
            self._load_results()

    def _on_process_change(self, _):
        pass

    def _on_monitor_counter_change(self, _):
        self._model.update_process_status_notification()

    def _on_load_results_click(self, _):
        self._load_results()

    def _load_results(self):
        self.children = [self.loading_message]
        self._render()
        self.rendered = True
        self._post_render()
        self.has_controls = False

    def _render_controls(self):
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

        self.children = [
            self.process_status_notification,
            ipw.HBox(
                children=[
                    self.load_results_button,
                    ipw.HTML("""
                        <div style="margin-left: 10px">
                            <b>Note:</b> Load time may vary depending on the size of the calculation
                        </div>
                    """),
                ]
            ),
        ]

        self.has_controls = True

    def _render(self):
        raise NotImplementedError()

    def _post_render(self):
        pass
