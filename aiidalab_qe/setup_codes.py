# -*- coding: utf-8 -*-
import os
from pathlib import Path
from subprocess import run
from threading import Thread

import ipywidgets as ipw
import traitlets
from aiida.common.exceptions import NotExistent
from aiida.orm import load_code
from filelock import FileLock, Timeout

from aiidalab_qe.widgets import ProgressBar

__all__ = [
    "QESetupWidget",
]

FN_LOCKFILE = Path.home().joinpath(".install-qe-on-localhost.lock")

QE_VERSION = "6.7"

CONDA_ENV_PREFIX = Path.home().joinpath(
    ".conda", "envs", f"quantum-espresso-{QE_VERSION}"
)

CODE_NAMES = ("pw", "projwfc", "dos")


def qe_installed():
    return CONDA_ENV_PREFIX.exists()


def install_qe():
    env = os.environ.copy()
    env["PATH"] = f"{env['PATH']}:{Path.home().joinpath('.local', 'bin')}"

    run(
        [
            "conda",
            "create",
            "--yes",
            "--prefix",
            str(CONDA_ENV_PREFIX),
            f"qe={QE_VERSION}",
        ],
        capture_output=True,
        check=True,
    )


def _code_is_setup(name):
    try:
        code = load_code(f"{name}-{QE_VERSION}@localhost")
        if (
            Path(code.get_remote_exec_path()).resolve()
            != CONDA_ENV_PREFIX.joinpath("bin", f"{name}.x").resolve()
        ):
            raise RuntimeError(f"Code {code} is already setup, but the paths differs!")
    except NotExistent:
        return False
    else:
        return True


def codes_are_setup():
    return all(_code_is_setup(code_name) for code_name in CODE_NAMES)


def _setup_code(code_name, computer_name="localhost"):
    try:
        load_code(f"{code_name}-{QE_VERSION}@localhost")
    except NotExistent:
        run(
            [
                "verdi",
                "code",
                "setup",
                "--non-interactive",
                "--label",
                f"{code_name}-{QE_VERSION}",
                "--description",
                f"{code_name}.x ({QE_VERSION}) setup by AiiDAlab.",
                "--input-plugin",
                f"quantumespresso.{code_name}",
                "--computer",
                computer_name,
                "--remote-abs-path",
                CONDA_ENV_PREFIX.joinpath("bin", f"{code_name}.x"),
            ],
            check=True,
            capture_output=True,
        )
    else:
        raise RuntimeError(f"Code {code_name} (v{QE_VERSION}) is already setup!")


def setup_codes():
    for code_name in CODE_NAMES:
        _setup_code(code_name)


class QESetupWidget(ProgressBar):

    installed = traitlets.Bool(allow_none=True).tag(readonly=True)
    busy = traitlets.Bool().tag(readonly=True)
    installing = traitlets.Bool().tag(readonly=True)
    error = traitlets.Unicode().tag(readonly=True)

    def __init__(self, prefix=None, hide_by_default=True, auto_start=True, **kwargs):
        self.prefix = prefix or f"QuantumESPRESSO (v{QE_VERSION}) @localhost: "
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
                    # Check whether the codes are present:
                    self.installed = codes_are_setup()
                    if not self.installed:
                        self.set_message("installing...")
                        # To setup our own codes, we install QE on the local
                        # host:
                        if not qe_installed():
                            self.set_message("Installing QE...")
                            self.value = ProgressBar.AnimationRate(0.05, max=0.7)
                            install_qe()
                        self.value = 0.7
                        # After installing QE, we install the corresponding
                        # AiiDA codes:
                        for i, code_name in enumerate(CODE_NAMES):
                            if not _code_is_setup(code_name):
                                self.set_message(
                                    f"Setting up AiiDA code ({code_name})..."
                                )
                                self.value = ProgressBar.AnimationRate(
                                    0.1, max=0.8 + i * 0.1
                                )
                                _setup_code(code_name)
                            self.value = 0.8 + i * 0.1
                        # After going through the installation procedure, we
                        # expect both our version of QE to be installed, as well
                        # as the codes to be setup.
                        self.installed = qe_installed() and codes_are_setup()

            except Timeout:
                # assume that the installation was triggered by a different
                # process
                self.set_message("installing...")
                self.value = self.AnimationRate(0.01)
                with FileLock(FN_LOCKFILE, timeout=120):
                    self.installed = codes_are_setup()

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
