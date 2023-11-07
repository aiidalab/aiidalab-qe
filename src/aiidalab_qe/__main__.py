"""For running the app from the command line used for post_install script.
"""

from pathlib import Path

import click
from aiida import load_profile

from aiidalab_qe.common.setup_codes import codes_are_setup
from aiidalab_qe.common.setup_codes import install as install_qe_codes

# The default profile name of AiiDAlab container.
_DEFAULT_PROFILE = "default"


@click.group()
def cli():
    pass


@cli.command()
@click.option("-f", "--force", is_flag=True)
@click.option("-p", "--profile", default=_DEFAULT_PROFILE)
def install_qe(force, profile):
    load_profile(profile)
    try:
        for msg in install_qe_codes(force=force):
            click.echo(msg)
        assert codes_are_setup()
        click.secho("Codes are setup!", fg="green")
    except Exception as error:
        raise click.ClickException(f"Failed to set up QE failed: {error}")


@cli.command()
@click.option("-p", "--profile", default=_DEFAULT_PROFILE, help="AiiDA profile name.")
@click.option(
    "source",
    "--source",
    default=None,
    help="The source folder to install from local.",
    type=click.Path(exists=True, path_type=Path, resolve_path=True),
)
def install_pseudos(profile, source):
    """Install pseudopotentials from a local folder if source is specified,
    otherwise download from remote repositories.
    """
    from aiidalab_qe.common.setup_pseudos import install

    load_profile(profile)

    try:
        for msg, _ in install(cwd=source, download_only=False):
            click.echo(msg)
        click.secho("Pseudopotentials are installed!", fg="green")
    except Exception as error:
        raise click.ClickException(f"Failed to set up pseudo potentials: {error}")


@cli.command()
@click.option(
    "dest",
    "--dest",
    default=None,
    help="The dest folder where to download the pseudos.",
    type=click.Path(exists=True, path_type=Path, resolve_path=True),
)
def download_pseudos(dest):
    from aiidalab_qe.common.setup_pseudos import EXPECTED_PSEUDOS, _install_pseudos

    try:
        for progress in _install_pseudos(
            EXPECTED_PSEUDOS, download_only=True, cwd=dest
        ):
            click.echo(progress)
        click.secho("Pseudopotentials are downloaded!", fg="green")

    except Exception as error:
        raise click.ClickException(f"Failed to download pseudo potentials: {error}")


if __name__ == "__main__":
    cli()
