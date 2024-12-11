from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl
from IPython.display import display

from aiidalab_qe.common.guide_manager import guide_manager
from aiidalab_qe.common.widgets import LoadingWidget


def without_triggering(toggle: str):
    """Decorator to prevent the other toggle from triggering its callback."""

    def decorator(func):
        def wrapper(self, change: dict):
            """Toggle off other button without triggering its callback."""
            view: AppWrapperView = self._view
            button: ipw.ToggleButton = getattr(view, toggle)
            callback = getattr(self, f"_on_{toggle}")
            button.unobserve(callback, "value")
            button.value = False
            func(self, change)
            button.observe(callback, "value")

        return wrapper

    return decorator


class AppWrapperContoller:
    """An MVC controller for `AppWrapper`."""

    def __init__(
        self,
        model: AppWrapperModel,
        view: AppWrapperView,
    ) -> None:
        """`AppWrapperController` constructor.

        Parameters
        ----------
        `model` : `AppWrapperModel`
            The associated model.
        `view` : `AppWrapperView`
            The associated view.
        """
        self._model = model
        self._view = view
        self._set_event_handlers()

    def enable_toggles(self) -> None:
        """Enable the toggle buttons."""
        self._view.guide_toggle.disabled = False
        self._view.about_toggle.disabled = False
        self._view.job_history_toggle.disabled = False

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
                    layout=ipw.Layout(align_items="baseline"),
                ),
            ]
            self._view.info_container.layout.display = "flex"
            self._view.job_history_toggle.value = False
        else:
            self._view.info_container.children = []
            self._view.info_container.layout.display = "none"

    @without_triggering("guide_toggle")
    def _on_about_toggle(self, change: dict):
        """Toggle the about section."""
        if change["new"]:
            self._view.info_container.children = [self._view.about]
            self._view.info_container.layout.display = "flex"
            self._view.job_history_toggle.value = False
        else:
            self._view.info_container.children = []
            self._view.info_container.layout.display = "none"

    def _on_job_history_toggle(self, change: dict):
        """Toggle the job list section."""
        if change["new"]:
            self._view.about_toggle.value = False
            self._view.guide_toggle.value = False
            self._old_view = self._view.main.children
            self._view.main.children = [LoadingWidget("Loading job history")]
            self._view.job_history.setup_table()
            self._view.main.children = [
                self._view.job_history.filters_layout,
                self._view.job_history.table,
            ]
        else:
            self._view.main.children = self._old_view

    def _on_guide_category_select(self, change: dict):
        self._view.guide_selection.options = guide_manager.get_guides(change["new"])
        self._update_active_guide()

    def _on_guide_select(self, _):
        self._update_active_guide()

    def _update_active_guide(self):
        """Sets the current active guide."""
        category = self._view.guide_category_selection.value
        guide = self._view.guide_selection.value
        active_guide = f"{category}/{guide}" if category != "none" else category
        guide_manager.active_guide = active_guide

    def _set_guide_category_options(self, _):
        """Fetch the available guides."""
        self._view.guide_category_selection.options = [
            "none",
            *guide_manager.get_guide_categories(),
        ]

    def _set_event_handlers(self) -> None:
        """Set up event handlers."""
        self._view.guide_toggle.observe(
            self._on_guide_toggle,
            "value",
        )
        self._view.about_toggle.observe(
            self._on_about_toggle,
            "value",
        )
        self._view.job_history_toggle.observe(
            self._on_job_history_toggle,
            "value",
        )
        self._view.guide_category_selection.observe(
            self._on_guide_category_select,
            "value",
        )
        self._view.guide_selection.observe(
            self._on_guide_select,
            "value",
        )
        self._view.on_displayed(self._set_guide_category_options)


class AppWrapperModel(tl.HasTraits):
    """An MVC model for `AppWrapper`."""

    def __init__(self):
        """`AppWrapperModel` constructor."""


class AppWrapperView(ipw.VBox):
    """An MVC view for `AppWrapper`."""

    def __init__(self) -> None:
        """`AppWrapperView` constructor."""

        ################# LAZY LOADING #################

        from datetime import datetime

        from importlib_resources import files
        from IPython.display import Image
        from jinja2 import Environment

        from aiidalab_qe.app.static import templates
        from aiidalab_qe.app.utils.search_jobs import QueryInterface
        from aiidalab_qe.common.infobox import InfoBox
        from aiidalab_qe.common.widgets import LoadingWidget
        from aiidalab_qe.version import __version__

        #################################################

        self.output = ipw.Output()

        logo_img = Image(
            filename="docs/source/_static/logo.png",
            width="700",
        )
        logo = ipw.Output()
        with logo:
            display(logo_img)
        logo.add_class("logo")

        subtitle = ipw.HTML("<h3 id='subtitle'>🎉 Happy computing 🎉</h3>")

        self.guide_toggle = ipw.ToggleButton(
            layout=ipw.Layout(width="auto"),
            button_style="",
            icon="book",
            value=False,
            description="Getting Started",
            tooltip="Learn how to use the app",
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

        self.job_history_toggle = ipw.ToggleButton(
            layout=ipw.Layout(width="auto"),
            button_style="",
            icon="list",
            value=False,
            description="Job History",
            tooltip="View all jobs run with this app",
            disabled=True,
        )

        info_toggles = ipw.HBox(
            children=[
                self.guide_toggle,
                self.about_toggle,
                self.job_history_toggle,
            ]
        )
        info_toggles.add_class("info-toggles")

        env = Environment()
        guide_template = files(templates).joinpath("guide.jinja").read_text()
        about_template = files(templates).joinpath("about.jinja").read_text()

        self.guide = ipw.HTML(env.from_string(guide_template).render())
        self.about = ipw.HTML(env.from_string(about_template).render())

        self.guide_category_selection = ipw.RadioButtons(
            options=["none"],
            description="Guides:",
            value="none",
            layout=ipw.Layout(width="max-content"),
        )
        self.guide_selection = ipw.RadioButtons(layout=ipw.Layout(margin="2px 20px"))

        self.job_history = QueryInterface()

        self.info_container = InfoBox()

        header = ipw.VBox(
            children=[
                logo,
                subtitle,
                info_toggles,
                self.info_container,
            ],
        )
        header.add_class("app-header")

        self.main = ipw.VBox(children=[LoadingWidget("Loading the app")])

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
                self.main,
                footer,
            ],
        )
