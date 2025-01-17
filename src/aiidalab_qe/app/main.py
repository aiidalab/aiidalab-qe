"""The main app that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

from pathlib import Path

from IPython.display import display

from aiidalab_qe.app.static import styles
from aiidalab_qe.app.wizard_app import WizardApp
from aiidalab_qe.app.wrapper import AppWrapperContoller, AppWrapperModel, AppWrapperView
from aiidalab_widgets_base.bug_report import (
    install_create_github_issue_exception_handler,
)
from aiidalab_widgets_base.utils.loaders import load_css


class QeApp:
    def __init__(self, process=None, qe_auto_setup=True, bug_report_url=None):
        """Initialize the AiiDAlab QE application with the necessary setup."""

        self.process = process
        self.qe_auto_setup = qe_auto_setup

        self._load_styles()

        # Initialize MVC components
        self.model = AppWrapperModel()
        self.view = AppWrapperView()
        display(self.view)

        # Set up bug report handling (if a URL is provided)
        if bug_report_url:
            install_create_github_issue_exception_handler(
                self.view.output,
                url=bug_report_url,
                labels=("bug", "automated-report"),
            )

        # setup UI controls
        self.controller = AppWrapperContoller(self.model, self.view)
        self.controller.enable_controls()

    def _load_styles(self):
        """Load CSS styles from the static directory."""
        load_css(css_path=Path(styles.__file__).parent)

    def load(self):
        """Initialize the WizardApp and integrate the app into the main view."""
        self.app = WizardApp(qe_auto_setup=self.qe_auto_setup)
        self.view.main.children = [self.app]
        # load a previous calculation if it is provided
        if self.process:
            self.app.process = self.process
