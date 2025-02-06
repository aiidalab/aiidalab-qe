from pathlib import Path
from threading import Thread

import ipywidgets as ipw
import traitlets

from ..setup.codes import QE_VERSION, install_and_setup
from .widgets import ProgressBar

__all__ = [
    "QESetupWidget",
]

FN_DO_NOT_SETUP = Path.cwd().joinpath(".do-not-setup-on-localhost")


class QESetupWidget(ipw.VBox):
    """The QE setup widget checks whether there are codes that match specific expected
    labels (e.g. "pw-7.4@localhost") and triggers both the installation of QE into a
    dedicated conda environment and the setup of the codes in case that they are not
    already configured.
    """

    installed = traitlets.Bool(allow_none=True).tag(readonly=True)
    busy = traitlets.Bool().tag(readonly=True)
    error = traitlets.Unicode().tag(readonly=True)

    def __init__(self, prefix=None, hide_by_default=True, auto_start=True, **kwargs):
        self.prefix = prefix or f"QuantumESPRESSO (v{QE_VERSION}) @localhost: "
        self.hide_by_default = hide_by_default

        self._progress_bar = ProgressBar(
            description=self.prefix,
            description_layout=ipw.Layout(min_width="300px"),
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
        )

        self._info_toggle_button = ipw.ToggleButton(
            icon="info-circle",
            disabled=True,
            layout=ipw.Layout(width="36px"),
        )
        self._info_toggle_button.observe(self._toggle_error_view, "value")

        self._reinstall_button = ipw.Button(
            icon="cogs",
            disabled=True,
            description="Install codes...",
            tooltip="Start another installation attempt.",
        )
        self._reinstall_button.on_click(self._trigger_reinstall)

        self._error_output = ipw.HTML()

        super().__init__(
            [
                ipw.HBox(
                    [self._progress_bar, self._info_toggle_button],
                    layout=ipw.Layout(width="auto"),
                ),
            ],
            **kwargs,
        )

        if auto_start:
            self.refresh()

    def set_message(self, msg):
        self._progress_bar.description = f"{self.prefix}{msg}"

    def _refresh_installed(self):
        try:
            self.set_trait("busy", True)

            for msg in install_and_setup():
                self.set_message(msg)

        except Exception as error:
            self.set_message("Failed to setup QE on localhost.")
            self.set_trait("error", str(error))
            FN_DO_NOT_SETUP.touch()
        else:
            self.set_trait("installed", True)
            self.set_message("OK")
        finally:
            self.set_trait("busy", False)

    def refresh(self):
        thread = Thread(target=self._refresh_installed)
        thread.start()

    @traitlets.default("installed")
    def _default_installed(self):
        return None

    @traitlets.default("busy")
    def _default_busy(self):
        return False

    @traitlets.default("failed")
    def _default_error(self):
        return ""

    @traitlets.observe("error")
    def _observe_error(self, change):
        with self.hold_trait_notifications():
            self._error_output.value = f"""
            <div class="alert alert-warning">
            <p>Failed to setup QE on localhost, due to error:</p>

            <p><code>{change["new"]}</code></p>

            <hr>
            <p>This means you have to setup QE manually to run it on this host.
            You can safely ignore this message if you do not plan on running
            QuantumESPRESSO calculations directly on the localhost. Alternatively
            you could try to make another installation attempt via the button
            below.</p>
            """
            self._info_toggle_button.disabled = not bool(change["new"])
            self._reinstall_button.disabled = not change["new"]
            if not change["new"]:
                self._info_toggle_button.value = False

    def _toggle_error_view(self, change):
        self.children = [self.children[0]] + (
            [self._error_output, self._reinstall_button] if change["new"] else []
        )

    @traitlets.observe("busy")
    @traitlets.observe("error")
    @traitlets.observe("installed")
    def _update(self, _change):
        with self.hold_trait_notifications():
            if self.hide_by_default:
                self.layout.visibility = (
                    "visible" if (self.busy or self.error) else "hidden"
                )

            if self.error or self.installed:
                self._progress_bar.value = 1.0
            elif self.busy:
                self._progress_bar.value = ProgressBar.AnimationRate(1.0)
            else:
                self._progress_bar.value = 0

            self._progress_bar.bar_style = (
                "info"
                if self.busy
                else (
                    "warning"
                    if self.error
                    else {True: "success", False: ""}.get(self.installed, "")
                )
            )

    def _trigger_reinstall(self, _=None):
        FN_DO_NOT_SETUP.unlink()
        self.refresh()
