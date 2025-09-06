"""The main app that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

from pathlib import Path

import ipywidgets as ipw
from IPython.display import display

from aiidalab_qe.app.static import styles
from aiidalab_widgets_base.bug_report import (
    install_create_github_issue_exception_handler,
)
from aiidalab_widgets_base.utils.loaders import load_css

from .app import AppController, AppModel, AppView

DEFAULT_BUG_REPORT_URL = "https://github.com/aiidalab/aiidalab-qe/issues/new"


class QeApp:
    def __init__(
        self,
        process=None,
        auto_setup=True,
        bug_report_url=DEFAULT_BUG_REPORT_URL,
        show_log=False,
    ):
        """Initialize the AiiDAlab QE application with the necessary setup."""

        self.process = process
        self.auto_setup = auto_setup
        self.log_widget = None

        self._load_styles()

        self.model = AppModel()
        self.view = AppView()
        display(self.view)

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
            display(
                ipw.VBox(
                    children=[
                        reset_button,
                        self.log_widget,
                    ],
                )
            )

        if bug_report_url:
            install_create_github_issue_exception_handler(
                self.log_widget if show_log else self.view.output,
                url=bug_report_url,
                labels=("bug", "automated-report"),
            )

        self.controller = AppController(self.model, self.view)
        self.controller.enable_toggles()

        if not self.model.validate_process(process):
            self.view.app_container.children = [
                ipw.HTML(f"""
                    <div class="alert alert-danger" style="text-align: center">
                        Process {process} does not exist
                        <br>
                        Please visit the <b>Calculation history</b> page to view
                        existing processes
                    </div>
                """)
            ]
        else:
            self.controller.load_app()

    def _load_styles(self):
        """Load CSS styles from the static directory."""
        load_css(css_path=Path(styles.__file__).parent)
