# -*- coding: utf-8 -*-
import os
from pathlib import Path
from subprocess import run
from threading import Thread

import ipywidgets as ipw
import traitlets
from aiida.orm import QueryBuilder
from aiida_pseudo.groups.family import PseudoPotentialFamily
from filelock import FileLock, Timeout

from aiidalab_qe.widgets import ProgressBar

EXPECTED_PSEUDOS = {
    "SSSP/1.1/PBE/efficiency",
    "SSSP/1.1/PBE/precision",
    "SSSP/1.1/PBEsol/efficiency",
    "SSSP/1.1/PBEsol/precision",
}


FN_LOCKFILE = Path.home().joinpath(".install-sssp.lock")


def pseudos_are_installed():
    qb = QueryBuilder()
    qb.append(PseudoPotentialFamily, filters={"label": {"like": "SSSP/%"}})
    labels = {group.label for group, in qb.iterall()}
    return labels == EXPECTED_PSEUDOS


def _install_pseudos():
    env = os.environ.copy()
    env["PATH"] = f"{env['PATH']}:{Path.home().joinpath('.local', 'bin')}"

    def run_(*args, **kwargs):
        return run(*args, env=env, capture_output=True, check=True, **kwargs)

    yield 0.1
    run_(["aiida-pseudo", "install", "sssp", "-x", "PBE", "-p", "efficiency"])
    yield 0.25
    run_(["aiida-pseudo", "install", "sssp", "-x", "PBE", "-p", "precision"])
    yield 0.5
    run_(["aiida-pseudo", "install", "sssp", "-x", "PBEsol", "-p", "efficiency"])
    yield 0.75
    run_(["aiida-pseudo", "install", "sssp", "-x", "PBEsol", "-p", "precision"])


class SSSPInstallWidget(ProgressBar):

    installed = traitlets.Bool(allow_none=True).tag(readonly=True)
    busy = traitlets.Bool().tag(readonly=True)
    installing = traitlets.Bool().tag(readonly=True)
    error = traitlets.Unicode().tag(readonly=True)

    def __init__(self, prefix=None, hide_by_default=True, auto_start=True, **kwargs):
        self.prefix = prefix or "SSSP pseudo potentials: "
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
        self.set_message("checking installation status...")
        try:
            self.set_trait("busy", True)
            try:
                with FileLock(FN_LOCKFILE, timeout=5):
                    self.installed = pseudos_are_installed()
                    if not self.installed:
                        self.set_message("installing...")
                        for progress in _install_pseudos():
                            self.value = progress
                        self.installed = pseudos_are_installed()

            except Timeout:
                # assume that the installation was triggered by a different process
                self.set_message("installing...")
                self.value = self.AnimationRate(0.01)
                with FileLock(FN_LOCKFILE, timeout=120):
                    self.installed = pseudos_are_installed()

            # Raise error in case that the installation was not successful
            # either in this process or a different one.
            if not self.installed:
                raise RuntimeError("Installation failed for unknown reasons.")

        except Exception as error:
            self.set_trait("error", str(error))
            self.set_message(str(error))
        else:
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
    def _update(self, change):
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
