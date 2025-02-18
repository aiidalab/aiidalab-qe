from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl
from IPython.display import display

from aiidalab_qe.app.wizard_app import WizardApp
from aiidalab_qe.common.guide_manager import guide_manager
from aiidalab_qe.common.setup_codes import QESetupWidget
from aiidalab_qe.common.setup_pseudos import PseudosInstallWidget
from aiidalab_qe.common.widgets import LinkButton, LoadingWidget
from aiidalab_widgets_base.bug_report import (
    install_create_github_issue_exception_handler,
)


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

    def setup(self):
        self._install_sssp()
        self._set_up_qe()

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

    def _on_installation_change(self, _):
        return
        if all([self._model.qe_installed, self._model.sssp_installed]):
            self.load_app()

    def _install_sssp(self):
        self.sssp_installation = PseudosInstallWidget(auto_start=False)
        self._view.main.children += (self.sssp_installation,)
        ipw.dlink(
            (self.sssp_installation, "installed"),
            (self._model, "sssp_installed"),
        )
        if self._model.qe_auto_setup:
            self.sssp_installation.refresh()

    def _set_up_qe(self):
        self.qe_setup = QESetupWidget(auto_start=False)
        self._view.main.children += (self.qe_setup,)
        ipw.dlink(
            (self.qe_setup, "installed"),
            (self._model, "qe_installed"),
        )
        if self._model.qe_auto_setup:
            self.qe_setup.refresh()

    def load_app(self):
        self._view.main.children = [LoadingWidget("Loading the app")]
        self.app = WizardApp(log_widget=self._view.log_widget)
        self._view.main.children = [self.app]
        self.app.process = self._model.process

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
        self._model.observe(
            self._on_installation_change,
            "sssp_installed",
        )
        self._model.observe(
            self._on_installation_change,
            "qe_installed",
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


class AppWrapperModel(tl.HasTraits):
    """An MVC model for `AppWrapper`."""

    process = tl.Union([tl.Unicode(), tl.Int()], allow_none=True)

    guide_category_options = tl.List(["none", *guide_manager.get_guide_categories()])
    selected_guide_category = tl.Unicode("none")
    guide_options = tl.List(tl.Unicode())
    selected_guide = tl.Unicode(None, allow_none=True)

    qe_installed = tl.Bool(allow_none=True)
    sssp_installed = tl.Bool(allow_none=True)

    qe_auto_setup = tl.Bool(True)
    show_log = tl.Bool(False)

    def update_active_guide(self, category, guide):
        """Sets the current active guide."""
        active_guide = f"{category}/{guide}" if category != "none" else category
        guide_manager.active_guide = active_guide


class AppWrapperView(ipw.VBox):
    """An MVC view for `AppWrapper`."""

    def __init__(self, show_log: bool = False, bug_report_url: str = "") -> None:
        """`AppWrapperView` constructor."""

        ################# LAZY LOADING #################

        from datetime import datetime

        from importlib_resources import files
        from IPython.display import Image
        from jinja2 import Environment

        from aiidalab_qe.app.static import images as images_folder
        from aiidalab_qe.app.static import templates
        from aiidalab_qe.common.infobox import InfoBox
        from aiidalab_qe.version import __version__

        #################################################

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

        self.guide_toggle = ipw.ToggleButton(
            layout=ipw.Layout(width="auto"),
            button_style="",
            icon="book",
            value=False,
            description="Getting started",
            tooltip="Learn how to use the app",
        )

        self.about_toggle = ipw.ToggleButton(
            layout=ipw.Layout(width="auto"),
            button_style="",
            icon="info",
            value=False,
            description="About",
            tooltip="Learn about the app",
        )

        self.calculation_history_link = LinkButton(
            description="Calculation history",
            link="./calculation_history.ipynb",
            icon="list",
        )

        self.setup_resources_link = LinkButton(
            description="Setup resources",
            link="../home/code_setup.ipynb",
            icon="database",
        )

        self.download_examples_link = LinkButton(
            description="Download examples",
            link="./examples.ipynb",
            icon="download",
        )

        self.new_workchain_link = LinkButton(
            description="New calculation",
            link="./qe.ipynb",
            icon="plus-circle",
        )

        self.controls = ipw.HBox(
            children=[
                self.guide_toggle,
                self.about_toggle,
                self.calculation_history_link,
                self.setup_resources_link,
                self.download_examples_link,
                self.new_workchain_link,
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
        self.guide_selection = ipw.RadioButtons(layout=ipw.Layout(margin="2px 20px"))

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

        self.main = ipw.VBox()

        current_year = datetime.now().year
        footer = ipw.HTML(f"""
            <footer>
                Copyright (c) {current_year} AiiDAlab team<br>
                Version: {__version__}
            </footer>
        """)

        self.log_container = ipw.VBox()

        if show_log:
            self.log_widget = ipw.Output(
                layout=ipw.Layout(
                    border="solid 1px lightgray",
                    margin="2px",
                    padding="5px",
                ),
            )
            reset_button = ipw.Button(
                description="Clear log",
                button_style="primary",
                icon="trash",
                layout=ipw.Layout(width="fit-content"),
            )
            reset_button.on_click(lambda _: self.log_widget.clear_output())
            self.log_container.children = [
                reset_button,
                self.log_widget,
            ]
        else:
            self.log_widget = None

        # Set up bug report handling (if a URL is provided)
        if bug_report_url:
            install_create_github_issue_exception_handler(
                self.log_widget if show_log else self.output,
                url=bug_report_url,
                labels=("bug", "automated-report"),
            )

        super().__init__(
            layout={},
            children=[
                self.output,
                header,
                self.main,
                footer,
                self.log_container,
            ],
        )
