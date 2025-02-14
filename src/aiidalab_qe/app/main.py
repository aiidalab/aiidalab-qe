"""The main app that shows the application in the Jupyter notebook.

Authors: AiiDAlab team
"""

from pathlib import Path

from IPython.display import display

from aiidalab_qe.app.static import styles
from aiidalab_qe.app.wrapper import AppWrapperContoller, AppWrapperModel, AppWrapperView
from aiidalab_widgets_base.utils.loaders import load_css

DEFAULT_BUG_REPORT_URL = "https://github.com/aiidalab/aiidalab-qe/issues/new"


class QeApp:
    def __init__(
        self,
        process=None,
        qe_auto_setup=True,
        bug_report_url=DEFAULT_BUG_REPORT_URL,
        show_log=False,
    ):
        """Initialize the AiiDAlab QE application with the necessary setup."""

        self._load_styles()

        # Initialize MVC components
        model = AppWrapperModel()
        model.process = process
        model.show_log = show_log
        model.qe_auto_setup = qe_auto_setup
        view = AppWrapperView(show_log, bug_report_url)
        _ = AppWrapperContoller(model, view)

        display(view)

    def _load_styles(self):
        """Load CSS styles from the static directory."""
        load_css(css_path=Path(styles.__file__).parent)
