from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from subprocess import run

from aiida_pseudo.groups.family import PseudoPotentialFamily
from filelock import FileLock, Timeout

from aiida.orm import QueryBuilder

SSSP_VERSION = "1.3"
PSEUDODOJO_VERSION = "0.4"

EXPECTED_PSEUDOS = {
    f"SSSP/{SSSP_VERSION}/PBE/efficiency",
    f"SSSP/{SSSP_VERSION}/PBE/precision",
    f"SSSP/{SSSP_VERSION}/PBEsol/efficiency",
    f"SSSP/{SSSP_VERSION}/PBEsol/precision",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/standard/upf",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/stringent/upf",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/stringent/upf",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/FR/standard/upf",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/FR/standard/upf",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/FR/stringent/upf",
    f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/FR/stringent/upf",
}


FN_LOCKFILE = Path.home().joinpath(".install-sssp.lock")


@dataclass
class PseudoFamily:
    """The dataclass to deal with pseudo family strings.

    Attributes:
    library: the library name of the pseudo family, e.g. SSSP or PseudoDojo.
    cmd_library_name: the sub command name used in aiida-pseudo command line.
    version: the version of the pseudo family, e.g. 1.2
    functional: the functional of the pseudo family, e.g. PBE, PBEsol.
    accuracy: the accuracy of the pseudo family, which is protocol in aiida-pseudo, e.g. efficiency, precision, standard, stringent.
    relativistic: the relativistic treatment of the pseudo family, e.g. SR, FR.
    file_type: the file type of the pseudo family, e.g. upf, psml, currently only used for PseudoDojo.
    """

    library: str
    version: str
    functional: str
    accuracy: str
    cmd_library_name: str = field(init=False)
    relativistic: str | None = None
    file_type: str | None = None

    def __post_init__(self):
        """Post init operations and checks."""
        if self.library == "SSSP":
            self.cmd_library_name = "sssp"
        elif self.library == "PseudoDojo":
            self.cmd_library_name = "pseudo-dojo"
        else:
            raise ValueError(f"Unknown pseudo library {self.library}")

    @classmethod
    def from_string(cls, pseudo_family_string: str) -> PseudoFamily:
        """Initialize from a pseudo family string."""
        # We support two pseudo families: SSSP and PseudoDojo
        # They are formatted as follows:
        # SSSP: SSSP/<version>/<functional>/<accuracy>
        # PseudoDojo: PseudoDojo/<version>/<functional>/<relativistic>/<accuracy>/<file_type>
        # where <relativistic> is either 'SR' or 'FR' and <file_type> is either 'upf' or 'psml'
        # Before we unify the format of family strings, the conditions below are necessary
        # to distinguish between the two families
        library = pseudo_family_string.split("/")[0]
        if library == "SSSP":
            version, functional, accuracy = pseudo_family_string.split("/")[1:]
            relativistic = None
            file_type = None
        elif library == "PseudoDojo":
            (
                version,
                functional,
                relativistic,
                accuracy,
                file_type,
            ) = pseudo_family_string.split("/")[1:]
        else:
            raise ValueError(
                f"Not able to parse valid library name from {pseudo_family_string}"
            )

        return cls(
            library=library,
            version=version,
            functional=functional,
            accuracy=accuracy,
            relativistic=relativistic,
            file_type=file_type,
        )


def pseudos_to_install() -> set[str]:
    """Query the database and return the list of pseudopotentials that are not installed."""
    qb = QueryBuilder()
    qb.append(
        PseudoPotentialFamily,
        filters={
            "or": [
                {"label": {"like": "SSSP/%"}},
                {"label": {"like": "PseudoDojo/%"}},
            ]
        },
        project="label",
    )
    labels = set(qb.all(flat=True))
    return EXPECTED_PSEUDOS - labels


def _construct_cmd(
    pseudo_family_string: str, download_only: bool = False, cwd: Path | None = None
) -> list:
    """Construct the command for installation of pseudopotentials.

    If ``cwd`` is not None, and ``download_only`` is True the, only download the
    pseudopotential files to the ``cwd`` folder.
    If ``download_only`` is False and ``cwd`` is not None, the the pseudos will be installed from the ``cwd`` where the pseudos are downloaded to.

    NOTE: download_only has nothing to do with cwd, it will not download the pseudos to cwd if cwd is specified.
    The control to download to cwd is in the ``_install_pseudos`` function below.
    """
    pseudo_family = PseudoFamily.from_string(pseudo_family_string)

    # the library used in command line is lowercase
    # e.g. SSSP -> sssp and PseudoDojo -> pseudo-dojo
    library = pseudo_family.cmd_library_name
    version = pseudo_family.version
    functional = pseudo_family.functional
    accuracy = pseudo_family.accuracy
    cmd = [
        "aiida-pseudo",
        "install",
        library,
        "--functional",
        functional,
        "--version",
        version,
        "-p",  # p for protocol which is the accuracy of the library
        accuracy,
    ]

    # extra arguments for PseudoDojo
    if library == "pseudo-dojo":
        relativistic = pseudo_family.relativistic
        file_type = pseudo_family.file_type
        cmd.extend(
            [
                "--relativistic",
                relativistic,
                "--pseudo-format",
                file_type,
            ]
        )

    if download_only:
        cmd.append("--download-only")

    # if cwd source folder specified, then install the pseudos from the folder
    # download file name is replace `/` with `_` of the pseudo family string with `.aiida_pseudo` extension
    if not download_only and cwd is not None:
        file_path = cwd / f"{pseudo_family_string.replace('/', '_')}.aiida_pseudo"
        if file_path.exists():
            cmd.extend(["--from-download", str(file_path)])

    return cmd


def run_cmd(cmd: list, env: dict | None = None, cwd: Path | None = None):
    """Run the command with specific env in the workdir specified."""
    run(cmd, env=env, cwd=cwd, capture_output=True, check=True)


def _install_pseudos(
    pseudo_families: set[str], download_only: bool = False, cwd: Path | None = None
) -> Iterable[float]:
    """Go through the list of pseudo families and install them."""
    env = os.environ.copy()
    env["PATH"] = f"{env['PATH']}:{Path.home() / '.local' / 'bin'}"

    mult = 1.0 / len(pseudo_families)
    yield mult * 0
    for i, pseudo_family in enumerate(pseudo_families):
        cmd = _construct_cmd(pseudo_family, download_only, cwd=cwd)

        run_cmd(cmd, env=env, cwd=cwd)

        yield mult * (i + 1)


def install(
    download_only: bool = False, cwd: Path | None = None
) -> Iterable[tuple[str, float]]:
    yield "Checking installation status...", 0.1
    try:
        with FileLock(FN_LOCKFILE, timeout=5):
            if len(pseudos := pseudos_to_install()) > 0:
                yield "Installing...", 0.1
                for progress in _install_pseudos(pseudos, download_only, cwd):
                    yield "Installing...", progress

    except Timeout:
        # Assume that the installation was triggered by a different process.
        from aiidalab_qe.common.widgets import ProgressBar

        yield (
            "Installation was already started elsewhere, waiting for it to finish...",
            ProgressBar.AnimationRate(1.0),
        )
        with FileLock(FN_LOCKFILE, timeout=300):
            if len(pseudos_to_install()) > 0:
                raise RuntimeError(
                    "Installation process did not finish in the expected time."
                ) from None
