from __future__ import annotations

from threading import Thread

import ipywidgets as ipw
import traitlets

from ..setup.pseudos import install, pseudos_to_install
from .widgets import ProgressBar


class PseudosInstallWidget(ProgressBar):
    """The SSSP installation status widget shows the installation status of the SSSP
    pseudo potentials and triggers the installation in case that they are not yet
    installed. The widget will remain in a "busy" state in case that the installation
    was already triggered elsewhere, e.g., by the start up scripts. The submission is
    blocked while the potentials are not yet installed.
    """

    installed = traitlets.Bool(allow_none=True).tag(readonly=True)
    busy = traitlets.Bool().tag(readonly=True)
    installing = traitlets.Bool().tag(readonly=True)
    error = traitlets.Unicode().tag(readonly=True)

    def __init__(self, prefix=None, hide_by_default=True, auto_start=True, **kwargs):
        self.prefix = prefix or "Pseudo potentials: "
        self.hide_by_default = hide_by_default

        super().__init__(
            description=self.prefix,
            description_layout=ipw.Layout(min_width="300px"),
            **kwargs,
        )

        if auto_start:
            self.refresh()

    def set_message(self, msg):
        self.description = f"{self.prefix}{msg}"

    def _refresh_installed(self):
        self.set_trait("busy", True)

        try:
            for msg, progress in install():
                self.set_message(msg)
                self.value = progress

        except Exception as error:
            self.set_trait("error", str(error))
            self.set_message(str(error))
        else:
            # If all the libraries are install by hands `pseudos_to_install()` will be empty list.
            self.set_trait("installed", not bool(pseudos_to_install()))
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
                self.value = 1.0

            self.bar_style = (
                "info"
                if self.busy
                else (
                    "danger"
                    if self.error
                    else {True: "success", False: "primary"}.get(self.installed, "")
                )
            )
