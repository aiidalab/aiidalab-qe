from __future__ import annotations

import ipywidgets as ipw
import traitlets


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
        self._view.job_list_toggle.disabled = False
        self._view.plugin_list_toggle.disabled = False

    @without_triggering("about_toggle")
    def _on_guide_toggle(self, change: dict):
        """Toggle the guide section."""
        self._view.info_container.children = [self._view.guide] if change["new"] else []
        self._view.info_container.layout.display = "flex" if change["new"] else "none"

    @without_triggering("guide_toggle")
    def _on_about_toggle(self, change: dict):
        """Toggle the about section."""
        self._view.info_container.children = [self._view.about] if change["new"] else []
        self._view.info_container.layout.display = "flex" if change["new"] else "none"

    @without_triggering("guide_toggle")
    def _on_job_list_toggle(self, change: dict):
        """Toggle the about section."""
        if change["new"]:
            self._view.job_list.setup_table()
            self._view.main.children = [
                self._view.job_list.filters_layout,
                self._view.job_list.table,
            ]
        else:
            self._view.main.children = [self._view.app]

    def _set_event_handlers(self) -> None:
        """Set up event handlers."""
        self._view.guide_toggle.observe(self._on_guide_toggle, "value")
        self._view.about_toggle.observe(self._on_about_toggle, "value")
        self._view.job_list_toggle.observe(self._on_job_list_toggle, "value")


class AppWrapperModel(traitlets.HasTraits):
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
        from IPython.display import Image, display
        from jinja2 import Environment

        from aiidalab_qe.app.static import templates
        from aiidalab_qe.app.utils.search_jobs import QueryInterface
        from aiidalab_qe.common.infobox import InfoBox
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
            button_style="",
            icon="question",
            value=False,
            description="Guide",
            tooltip="Learn how to use the app",
            disabled=True,
        )

        self.about_toggle = ipw.ToggleButton(
            button_style="",
            icon="info",
            value=False,
            description="About",
            tooltip="Learn about the app",
            disabled=True,
        )

        self.job_list_toggle = ipw.ToggleButton(
            button_style="",
            icon="info",
            value=False,
            description="Job List",
            tooltip="Learn about the app",
            disabled=True,
        )

        self.plugin_list_toggle = ipw.ToggleButton(
            button_style="",
            icon="info",
            value=False,
            description="Plugin List",
            tooltip="Learn about the app",
            disabled=True,
        )

        info_toggles = ipw.HBox(
            children=[
                self.guide_toggle,
                self.about_toggle,
                self.job_list_toggle,
                self.plugin_list_toggle,
            ]
        )
        info_toggles.add_class("info-toggles")

        env = Environment()
        guide_template = files(templates).joinpath("guide.jinja").read_text()
        about_template = files(templates).joinpath("about.jinja").read_text()

        self.guide = ipw.HTML(env.from_string(guide_template).render())
        self.about = ipw.HTML(env.from_string(about_template).render())

        self.info_container = InfoBox()
        self.job_list = QueryInterface()

        header = ipw.VBox(
            children=[
                logo,
                subtitle,
                info_toggles,
                self.info_container,
            ],
        )
        header.add_class("app-header")

        loading = ipw.HTML("""
            <div id="loading">
                Loading the app <i class="fa fa-spinner fa-spin"></i>
            </div>
        """)

        self.main = ipw.VBox(children=[loading])

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
