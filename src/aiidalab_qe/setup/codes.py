import json
import subprocess
from pathlib import Path
from shutil import which

from filelock import FileLock, Timeout

from aiida.common.exceptions import NotExistent
from aiida.orm import load_code

FN_INSTALL_LOCKFILE = Path.home().joinpath(".install-qe-on-localhost.lock")
FN_SETUP_LOCKFILE = Path.home().joinpath(".setup-qe-on-localhost.lock")
FN_DO_NOT_SETUP = Path.cwd().joinpath(".do-not-setup-on-localhost")

QE_VERSION = "7.4"


def get_qe_env():
    # QE is already pre-installed in the QE image
    path = Path(f"/opt/conda/envs/quantum-espresso-{QE_VERSION}")
    if path.exists():
        return path
    else:
        return Path.home().joinpath(".conda", "envs", f"quantum-espresso-{QE_VERSION}")


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
    """Check if Quantum Espresso (QE) is installed in the specified conda environment.

    Returns:
        bool: True if the environment exists and QE is installed; False otherwise.
    """
    try:
        # Verify if the specified conda environment exists
        env_exist = get_qe_env().exists()

        if not env_exist:
            return False

        # Run the conda list command to check for the QE package
        proc = subprocess.run(
            [
                "conda",
                "list",
                "-n",
                f"{get_qe_env().name}",
                "--json",
                "--full-name",
                "qe",
            ],
            check=True,
            capture_output=True,
        )

        # Load and interpret the JSON output
        info = json.loads(proc.stdout.decode())

        # Check if 'qe' is listed in the environment
        for package in info:
            if package.get("name") == "qe":
                return True
        return False  # noqa: TRY300
    except Exception as error:
        raise RuntimeError(
            "Failed to check if Quantum Espresso is installed."
        ) from error


def install_qe():
    subprocess.run(
        [
            "conda",
            "create",
            "--yes",
            "--override-channels",
            "--channel",
            "conda-forge",
            "--prefix",
            str(get_qe_env()),
            f"qe={QE_VERSION}",
        ],
        capture_output=True,
        check=True,
    )


def _code_is_setup(name, computer):
    try:
        load_code(f"{name}-{QE_VERSION}@{computer}")
    except NotExistent:
        return False
    else:
        return True


def codes_are_setup(computer):
    return all(_code_is_setup(code_name, computer) for code_name in CODE_NAMES)


def _generate_header_to_setup_code():
    """Generate the header string to setup a code for a given computer."""
    header_code = """
from aiida.orm.nodes.data.code.installed import InstalledCode
from aiida.orm import load_computer
from aiida import load_profile
load_profile()

"""
    return header_code


def _generate_string_to_setup_code(code_name, computer):
    """Generate the Python string to setup an AiiDA code for a given computer.

    Tries to load an existing code and if not existent,
    generates Python code to create and store a new code setup."""
    try:
        load_code(f"{code_name}-{QE_VERSION}@{computer}")
    except NotExistent:
        label = f"{code_name}-{QE_VERSION}"
        description = f"{code_name}.x ({QE_VERSION}) setup by AiiDAlab."
        filepath_executable = get_qe_env().joinpath("bin", f"{code_name}.x")
        default_calc_job_plugin = f"quantumespresso.{code_name}"
        prepend_text = f'eval "$(conda shell.posix hook)"\\nconda activate {get_qe_env()}\\nexport OMP_NUM_THREADS=1'
        python_code = """
computer = load_computer('{}')
code = InstalledCode(computer=computer,
                    label='{}',
                    description='{}',
                    filepath_executable='{}',
                    default_calc_job_plugin='{}',
                    prepend_text='{}'
                    )

code.store()
""".format(  # noqa: UP032
            computer,
            label,
            description,
            filepath_executable,
            default_calc_job_plugin,
            prepend_text,
        )
        return python_code
    else:
        # the code already exists
        return ""


def setup_codes(computer):
    python_code = _generate_header_to_setup_code()
    for code_name in CODE_NAMES:
        python_code += _generate_string_to_setup_code(code_name, computer)
    try:
        subprocess.run(["python", "-c", python_code], capture_output=True, check=True)
    except subprocess.CalledProcessError as err:
        raise RuntimeError(
            f"Failed to setup codes, exit_code={err.returncode}, {err.stderr}"
        ) from None


def install_and_setup(computer="localhost", force=False):
    """Install Quantum ESPRESSO and the corresponding AiiDA codes.

    Args:
        force: Ignore previously failed attempts and install anyways.
        computer: computer label in AiiDA where the code is setup for
    """
    # Check for "do not install file" and skip actual check. The purpose of
    # this file is to not re-try this process on every app start in case that
    # there are issues.
    # XXX: use filelock to control `FN_DO_NOT_SETUP` as well
    if not force and FN_DO_NOT_SETUP.exists():
        raise RuntimeError("Installation failed in previous attempt.")

    yield from _install()
    yield from _setup(computer)


def _install():
    """Install Quantum ESPRESSO."""
    yield "Checking installation status..."

    conda_installed = which("conda")
    try:
        with FileLock(FN_INSTALL_LOCKFILE, timeout=5):
            if not conda_installed:
                raise RuntimeError(
                    "Unable to automatically install Quantum ESPRESSO, conda "
                    "is not available."
                )

            if qe_installed():
                return

            # Install Quantum ESPRESSO.
            yield "Installing QE..."
            try:
                install_qe()
            except subprocess.CalledProcessError as error:
                raise RuntimeError(
                    f"Failed to create conda environment: {error}"
                ) from None

    except Timeout:
        # Assume that the installation was triggered by a different process.
        yield "Installation was already started, waiting for it to finish..."
        with FileLock(FN_INSTALL_LOCKFILE, timeout=120):
            if not qe_installed():
                raise RuntimeError(
                    "Installation process did not finish in the expected time."
                ) from None


def _setup(computer):
    """Setup the corresponding AiiDA codes after QE installation."""
    yield "Checking setup status..."

    try:
        with FileLock(FN_SETUP_LOCKFILE, timeout=5):
            # We assume that if the codes are already setup, everything is in
            # order. Only if they are not present, should we take action,
            # however we only do so if the environment has a conda binary
            # present (`which conda`). If that is not the case then we assume
            # that this is a custom user environment in which case we also take
            # no further action.
            if codes_are_setup(computer=computer):
                return  # Already setup

            # After installing QE, we install the corresponding
            # AiiDA codes:
            python_code = _generate_header_to_setup_code()
            for code_name in CODE_NAMES:
                if not _code_is_setup(code_name, computer=computer):
                    yield f"Preparing setup script for ({code_name}) on ({computer})..."
                    code_string = _generate_string_to_setup_code(code_name, computer)
                    python_code += code_string
            try:
                yield "Setting up all codes..."
                subprocess.run(
                    ["python", "-c", python_code], capture_output=True, check=True
                )
            except subprocess.CalledProcessError as err:
                raise RuntimeError(
                    f"Failed to setup codes, exit_code={err.returncode}, {err.stderr}"
                ) from None

    except Timeout:
        # Assume that the installation was triggered by a different process.
        yield "Installation was already started, waiting for it to finish..."
        with FileLock(FN_SETUP_LOCKFILE, timeout=120):
            if not codes_are_setup(computer=computer):
                raise RuntimeError(
                    "Installation process did not finish in the expected time."
                ) from None
