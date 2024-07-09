from __future__ import annotations

from datetime import datetime

import ipywidgets as ipw
import traitlets
from aiidalab_widgets_base.infobox import InfoBox
from importlib_resources import files
from IPython.display import Image, Javascript, display
from jinja2 import Environment

from aiidalab_qe.app import static
from aiidalab_qe.version import __version__


def without_triggering(toggle: str):
    """Decorator to prevent the other toggle from triggering its callback."""

    def decorator(func):
        def wrapper(self, change: dict):
            """Toggle off other button without triggering its callback."""
            view: AppWrapperView = getattr(self, "_view")
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

    @without_triggering("about_toggle")
    def _on_guide_toggle(self, change: dict):
        """Toggle the guide section."""
        self._view.info_container.children = (
            [
                self._view.guide,
                self._view.guide_selection,
            ]
            if change["new"]
            else []
        )
        self._view.info_container.layout.display = "flex" if change["new"] else "none"

    @without_triggering("guide_toggle")
    def _on_about_toggle(self, change: dict):
        """Toggle the about section."""
        self._view.info_container.children = [self._view.about] if change["new"] else []
        self._view.info_container.layout.display = "flex" if change["new"] else "none"

    def _on_guide_select(self, change: dict):
        """Toggle the guide section."""
        display(
            Javascript(f"""
                document.querySelectorAll('.{change["old"]}').forEach((guide) => {'{'}
                    guide.classList.remove('show');
                {'}'});
            """)
        )
        if (guide_class := change["new"]) != "none":
            display(
                Javascript(f"""
                    document.querySelectorAll('.{guide_class}').forEach((guide) => {'{'}
                        guide.classList.add('show');
                    {'}'});
                """)
            )

    def _on_close_first_time_info(self, _=None):
        """Close the first time info box."""
        self._view.first_time_users_infobox.layout.display = "none"
        with open("first-time-user", "w") as file:
            file.write("existing user")

    def _set_event_handlers(self) -> None:
        """Set up event handlers."""
        self._view.guide_toggle.observe(self._on_guide_toggle, "value")
        self._view.about_toggle.observe(self._on_about_toggle, "value")
        self._view.guide_selection.observe(self._on_guide_select, "value")
        self._view.close_first_time_info_button.on_click(self._on_close_first_time_info)


class AppWrapperModel(traitlets.HasTraits):
    """An MVC model for `AppWrapper`."""

    def __init__(self):
        """`AppWrapperModel` constructor."""


class AppWrapperView(ipw.VBox):
    """An MVC view for `AppWrapper`."""

    def __init__(self) -> None:
        """`AppWrapperView` constructor."""
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

        self.close_first_time_info_button = ipw.Button(
            icon="times",
            tooltip="Close",
        )

        self.first_time_users_infobox = InfoBox(
            title="First time users",
            children=[
                self.close_first_time_info_button,
                ipw.HTML("""
                    <p>
                        If you are <strong>new to the Quantum ESPRESSO app</strong>,
                        click on the guide button below to learn the basics of the
                        app and/or follow along to the tutorials.
                    </p>
                """),
            ],
            **{"custom-css": "first-time-users-infobox"},
        )

        self._check_if_first_time_user()

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

        info_toggles = ipw.HBox(
            children=[
                self.guide_toggle,
                self.about_toggle,
            ]
        )
        info_toggles.add_class("info-toggles")

        env = Environment()
        guide_template = files(static).joinpath("guide.jinja").read_text()
        about_template = files(static).joinpath("about.jinja").read_text()

        self.guide = ipw.HTML(env.from_string(guide_template).render())
        self.about = ipw.HTML(env.from_string(about_template).render())

        self.guide_selection = ipw.RadioButtons(
            options=[
                "none",
                "relaxation",
                "bands",
            ],
            description="Guides",
            value="none",
        )

        self.info_container = InfoBox()

        header = ipw.VBox(
            children=[
                logo,
                subtitle,
                self.first_time_users_infobox,
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

    def _check_if_first_time_user(self):
        """Check if the user is a first time user."""
        try:
            with open("first-time-user") as file:
                first_time_user = file.read().find("existing user") == -1
        except FileNotFoundError:
            first_time_user = True

        self.first_time_users_infobox.layout.display = (
            "flex" if first_time_user else "none"
        )
