from __future__ import annotations

import typing as t
from datetime import datetime

import ipywidgets as ipw
import traitlets as tl
from importlib_resources import files
from IPython.display import Image, display
from jinja2 import Environment

from aiida import orm
from aiida.orm.utils.serialize import deserialize_unsafe
from aiidalab_qe.app.static import images as images_folder
from aiidalab_qe.app.static import templates
from aiidalab_qe.common.guide_manager import guide_manager
from aiidalab_qe.common.infobox import InAppGuide, InfoBox
from aiidalab_qe.common.widgets import LinkButton
from aiidalab_qe.version import __version__
from aiidalab_widgets_base import LoadingWidget

from .wizard import Wizard, WizardModel


def without_triggering(toggle: str):
    """Decorator to prevent the other toggle from triggering its callback."""

    def decorator(func):
        def wrapper(self, change: dict):
            """Toggle off other button without triggering its callback."""
            view: AppView = self._view
            button: ipw.ToggleButton = getattr(view, toggle)
            callback = getattr(self, f"_on_{toggle}")
            button.unobserve(callback, "value")
            button.value = False
            func(self, change)
            button.observe(callback, "value")

        return wrapper

    return decorator


class AppController:
    """An MVC controller for the app."""

    def __init__(self, model: AppModel, view: AppView) -> None:
        """`AppController` constructor.

        Parameters
        ----------
        `model` : `AppModel`
            The associated model.
        `view` : `AppView`
            The associated view.
        """
        self._model = model
        self._view = view
        self._wizard_model = WizardModel()
        self._set_event_handlers()

    def enable_toggles(self) -> None:
        """Enable the toggle buttons at the top of the app."""
        toggle: ipw.ToggleButton
        for toggle in self._view.toggles.children:
            toggle.disabled = False

    def load_wizard(
        self,
        auto_setup: bool = True,
        log_widget: ipw.Output | None = None,
    ) -> None:
        """Load and initialize the wizard."""
        self.wizard = Wizard(self._wizard_model, auto_setup, log_widget)
        self._view.app_container.children = [self.wizard]
        state = {"process_uuid": self._model.process_uuid}
        if self._model.process_uuid:
            state |= self._model.get_state_from_process()
        self._wizard_model.state = state
        self._model.loaded = True

    @without_triggering("about_toggle")
    def _on_guide_toggle(self, change: dict):
        """Toggle the guide section."""
        if change["new"]:
            self._view.info_container.children = [
                self._view.guide,
                ipw.HBox(
                    children=[
                        self._view.guide_category_selection,
                        self._view.guide_selection,
                    ],
                    layout=ipw.Layout(align_items="baseline", grid_gap="10px"),
                ),
            ]
            self._view.info_container.layout.display = "flex"
        else:
            self._view.info_container.children = []
            self._view.info_container.layout.display = "none"

    @without_triggering("guide_toggle")
    def _on_about_toggle(self, change: dict):
        """Toggle the about section."""
        if change["new"]:
            self._view.info_container.children = [self._view.about]
            self._view.info_container.layout.display = "flex"
        else:
            self._view.info_container.children = []
            self._view.info_container.layout.display = "none"

    def _on_guide_category_selection_change(self, change):
        self._model.guide_options = guide_manager.get_guides(change["new"])

    def _on_guide_selection_change(self, _):
        category = self._view.guide_category_selection.value
        guide = self._view.guide_selection.value
        self._model.update_active_guide(category, guide)

    def _set_event_handlers(self) -> None:
        """Set up event handlers."""
        self._model.observe(
            self._on_guide_category_selection_change,
            "selected_guide_category",
        )
        self._model.observe(
            self._on_guide_selection_change,
            [
                "selected_guide_category",
                "selected_guide",
            ],
        )

        self._view.guide_toggle.observe(
            self._on_guide_toggle,
            "value",
        )
        self._view.about_toggle.observe(
            self._on_about_toggle,
            "value",
        )

        ipw.dlink(
            (self._model, "guide_category_options"),
            (self._view.guide_category_selection, "options"),
        )
        ipw.link(
            (self._model, "selected_guide_category"),
            (self._view.guide_category_selection, "value"),
        )
        ipw.dlink(
            (self._model, "guide_options"),
            (self._view.guide_selection, "options"),
        )
        ipw.link(
            (self._model, "selected_guide"),
            (self._view.guide_selection, "value"),
        )
        ipw.dlink(
            (self._wizard_model, "loading"),
            (self._view.wizard_loading_message.layout, "display"),
            lambda loading: "flex" if loading else "none",
        )


class AppModel(tl.HasTraits):
    """An MVC model for the app."""

    guide_category_options = tl.List(
        ["No guides", *guide_manager.get_guide_categories()]
    )
    selected_guide_category = tl.Unicode("No guides")
    guide_options = tl.List(tl.Unicode())
    selected_guide = tl.Unicode(None, allow_none=True)

    process_uuid: str | None = None
    loaded = False

    @property
    def process(self) -> orm.WorkChainNode:
        return t.cast(orm.WorkChainNode, orm.load_node(self.process_uuid))

    def validate_process(self, process_identifier: str | None = None) -> bool:
        """Validate the process identifier."""
        if process_identifier:
            try:
                process = orm.load_node(process_identifier)
                assert isinstance(process, orm.WorkChainNode)
                self.process_uuid = process.uuid
            except Exception:
                return False
        return True

    def update_active_guide(self, category, guide):
        """Sets the current active guide."""
        active_guide = f"{category}/{guide}" if category != "No guides" else category
        guide_manager.active_guide = active_guide

    def get_state_from_process(self) -> dict:
        if not self.process_uuid:
            return {}
        parameters: dict = self.process.base.extras.get("ui_parameters", {})
        if parameters and isinstance(parameters, str):
            parameters = deserialize_unsafe(parameters)
        codes = parameters.pop("codes", {})

        # BACKWARDS COMPATIBILITY
        # We used to store the codes under "resources"
        if "resources" in parameters:
            resources = parameters["resources"]
            codes |= {key: {"code": value} for key, value in codes.items()}
            codes["pw"]["nodes"] = resources["num_machines"]
            codes["pw"]["cpus"] = resources["num_mpiprocs_per_machine"]
            codes["pw"]["parallelization"] = {"npool": resources["npools"]}
        # END BACKWARDS COMPATIBILITY

        return {
            "structure_state": {"uuid": self.process.inputs.structure.uuid},
            "configuration_state": parameters,
            "resources_state": codes,
        }


class AppView(ipw.VBox):
    """An MVC view for the app."""

    def __init__(self) -> None:
        """`AppView` constructor."""

        self.output = ipw.Output()

        logo_img = Image(
            filename=files(images_folder) / "logo.png",
            width="700",
        )
        logo = ipw.Output()
        with logo:
            display(logo_img)
        logo.add_class("logo")

        subtitle = ipw.HTML("<h3 id='subtitle'>🎉 Happy computing 🎉</h3>")

        self.calculation_history_link = LinkButton(
            description="Calculation history",
            link="./calculation_history.ipynb",
            icon="list",
            tooltip="View a list of previous calculations",
            style_="background-color: var(--color-aiida-orange)",
        )

        self.setup_resources_link = LinkButton(
            description="Setup resources",
            link="../home/code_setup.ipynb",
            icon="database",
            tooltip="Setup computational resources for your calculations",
            style_="background-color: var(--color-aiida-blue)",
        )

        self.new_workchain_link = LinkButton(
            description="New calculation",
            link="./qe.ipynb",
            icon="plus-circle",
            tooltip="Open a new calculation in a separate tab",
            style_="background-color: var(--color-aiida-green)",
        )

        self.external_links = ipw.HBox(
            children=[
                self.calculation_history_link,
                self.setup_resources_link,
                self.new_workchain_link,
            ],
        )

        self.guide_toggle = ipw.ToggleButton(
            layout=ipw.Layout(width="auto"),
            button_style="",
            icon="book",
            value=False,
            description="Tutorials",
            tooltip="Learn how to use the app through dedicated in-app guides",
            disabled=True,
        )

        self.about_toggle = ipw.ToggleButton(
            layout=ipw.Layout(width="auto"),
            button_style="",
            icon="info",
            value=False,
            description="About",
            tooltip="Learn about the app",
            disabled=True,
        )

        self.external_links.add_class("app-external-links")

        self.toggles = ipw.HBox(
            children=[
                self.guide_toggle,
                self.about_toggle,
            ],
        )
        self.toggles.add_class("app-toggles")

        self.controls = ipw.Box(
            children=[
                self.external_links,
                self.toggles,
            ],
        )
        self.controls.add_class("app-controls")

        env = Environment()
        guide_template = files(templates).joinpath("guide.jinja").read_text()
        about_template = files(templates).joinpath("about.jinja").read_text()

        self.guide = ipw.HTML(env.from_string(guide_template).render())
        self.about = ipw.HTML(env.from_string(about_template).render())

        self.guide_category_selection = ipw.RadioButtons(
            description="Guides:",
            layout=ipw.Layout(width="max-content"),
        )
        self.guide_selection = ipw.RadioButtons(layout=ipw.Layout(width="max-content"))

        self.info_container = InfoBox(layout=ipw.Layout(margin="14px 2px 0"))

        header = ipw.VBox(
            children=[
                logo,
                subtitle,
                self.controls,
                self.info_container,
            ],
        )
        header.add_class("app-header")

        self.wizard_loading_message = LoadingWidget(
            message="Populating the wizard",
            layout={"display": "none"},
        )

        self.app_container = ipw.VBox(children=[LoadingWidget("Loading the app")])

        current_year = datetime.now().year
        footer = ipw.HTML(f"""
            <footer>
                Copyright (c) {current_year} AiiDAlab team<br>
                Version: {__version__}
            </footer>
        """)

        super().__init__(
            layout={},
            children=[
                self.output,
                header,
                InAppGuide(identifier="guide-header"),
                self.wizard_loading_message,
                self.app_container,
                InAppGuide(identifier="post-guide"),
                footer,
            ],
        )
