# -*- coding: utf-8 -*-
from pathlib import Path
from shutil import which
from subprocess import CalledProcessError, run
from threading import Thread

import ipywidgets as ipw
import traitlets
from aiida.common.exceptions import NotExistent
from aiida.orm import load_code
from filelock import FileLock, Timeout

from aiidalab_qe.common.widgets import ProgressBar

__all__ = [
    "QESetupWidget",
]

FN_LOCKFILE = Path.home().joinpath(".install-qe-on-localhost.lock")
FN_DO_NOT_SETUP = Path.cwd().joinpath(".do-not-setup-on-localhost")

QE_VERSION = "7.2"

CONDA_ENV_PREFIX = Path.home().joinpath(
    ".conda", "envs", f"quantum-espresso-{QE_VERSION}"
)

# Add all QE codes with the calcjob entry point in the aiida-quantumespresso.
CODE_NAMES = (
    "pw",
    "projwfc",
    "dos",
    "cp",
    "epw",
    "matdyn",
    "neb",
    "open_grid",
    "ph",
    "pp",
    "pw2gw",
    "pw2wannier90",
    "q2r",
    "xspectra",
    "hp",
)


def qe_installed():
    return CONDA_ENV_PREFIX.exists()


def install_qe():
    run(
        [
            "conda",
            "create",
            "--yes",
            "--override-channels",
            "--channel",
            "conda-forge",
            "--prefix",
            str(CONDA_ENV_PREFIX),
            f"qe={QE_VERSION}",
        ],
        capture_output=True,
        check=True,
    )


def _code_is_setup(name):
    try:
        load_code(f"{name}-{QE_VERSION}@localhost")
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
                "create",
                "core.code.installed",
                "--non-interactive",
                "--label",
                f"{code_name}-{QE_VERSION}",
                "--description",
                f"{code_name}.x ({QE_VERSION}) setup by AiiDAlab.",
                "--default-calc-job-plugin",
                f"quantumespresso.{code_name}",
                "--computer",
                computer_name,
                "--prepend-text",
                f'eval "$(conda shell.posix hook)"\nconda activate {CONDA_ENV_PREFIX}\nexport OMP_NUM_THREADS=1',
                "--filepath-executable",
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


def install(force=False):
    """Install Quantum ESPRESSO and the corresponding AiiDA codes.

    Args:
        force: Ignore previously failed attempts and install anyways.
    """
    # Check for "do not install file" and skip actual check. The purpose of
    # this file is to not re-try this process on every app start in case that
    # there are issues.
    if not force and FN_DO_NOT_SETUP.exists():
        raise RuntimeError("Installation failed in previous attempt.")

    yield "Checking installation status..."

    conda_installed = which("conda")
    try:
        with FileLock(FN_LOCKFILE, timeout=5):
            # We assume that if the codes are already setup, everything is in
            # order. Only if they are not present, should we take action,
            # however we only do so if the environment has a conda binary
            # present (`which conda`). If that is not the case then we assume
            # that this is a custom user environment in which case we also take
            # no further action.
            if codes_are_setup():
                return  # Already setup

            if not conda_installed:
                raise RuntimeError(
                    "Unable to automatically install Quantum ESPRESSO, conda "
                    "is not available."
                )

            if not qe_installed():
                # First, install Quantum ESPRESSO.
                yield "Installing QE..."
                try:
                    install_qe()
                except CalledProcessError as error:
                    raise RuntimeError(f"Failed to create conda environment: {error}")

            # After installing QE, we install the corresponding
            # AiiDA codes:
            for code_name in CODE_NAMES:
                if not _code_is_setup(code_name):
                    yield f"Setting up AiiDA code ({code_name})..."
                    _setup_code(code_name)

    except Timeout:
        # Assume that the installation was triggered by a different process.
        yield "Installation was already started, waiting for it to finish..."
        with FileLock(FN_LOCKFILE, timeout=120):
            if not codes_are_setup():
                raise RuntimeError(
                    "Installation process did not finish in the expected time."
                )


class QESetupWidget(ipw.VBox):
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

            for msg in install():
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
    def _update(self, change):
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
