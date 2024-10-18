"""For running the app from the command line used for post_install script."""

import sys
from pathlib import Path

import click

# The default profile name of AiiDAlab container.
_DEFAULT_PROFILE = "default"


@click.group()
def cli():
    pass


@cli.command()
@click.option("-f", "--force", is_flag=True)
@click.option("--computer")
@click.option("-p", "--profile", default=_DEFAULT_PROFILE)
def install_qe(force, profile, computer):
    from aiida import load_profile
    from aiidalab_qe.setup.codes import codes_are_setup, install_and_setup

    load_profile(profile)
    try:
        for msg in install_and_setup(computer=computer, force=force):
            click.echo(msg)
        assert codes_are_setup(computer=computer)
        click.secho("Codes are setup!", fg="green")
    except Exception as error:
        raise click.ClickException(f"Failed to set up QE: {error}") from error


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
    from aiida import load_profile
    from aiidalab_qe.setup.pseudos import install

    load_profile(profile)

    try:
        for msg, _ in install(cwd=source, download_only=False):
            click.echo(msg)
        click.secho("Pseudopotentials are installed!", fg="green")
    except Exception as error:
        raise click.ClickException(
            f"Failed to set up pseudo potentials: {error}"
        ) from error


@cli.command()
@click.option(
    "dest",
    "--dest",
    default=None,
    help="The dest folder where to download the pseudos.",
    type=click.Path(exists=True, path_type=Path, resolve_path=True),
)
def download_pseudos(dest):
    from aiidalab_qe.setup.pseudos import EXPECTED_PSEUDOS, _install_pseudos

    try:
        for progress in _install_pseudos(
            EXPECTED_PSEUDOS, download_only=True, cwd=dest
        ):
            click.echo(progress)
        click.secho("Pseudopotentials are downloaded!", fg="green")

    except Exception as error:
        raise click.ClickException(
            f"Failed to download pseudo potentials: {error}"
        ) from error


@cli.command()
@click.argument(
    "plugin_name",
    default="aiidalab_qe",
)
@click.option("-p", "--profile", default=_DEFAULT_PROFILE)
def test_plugin(plugin_name, profile):
    from aiida import load_profile
    from aiidalab_qe.app.utils import test_plugin_functionality

    load_profile(profile)

    try:
        success, message = test_plugin_functionality(plugin_name)
        if success:
            click.secho("Plugin is loaded successfully!", fg="green")
        else:
            click.secho(f"Failed to load plugin: {message}", fg="red", err=True)
            sys.exit(1)  # Exit with status 1 to indicate failure
    except Exception as error:
        click.secho(f"Failed to load plugin: {error}", fg="red", err=True)
        sys.exit(1)  # Exit with status 1 to indicate failure


if __name__ == "__main__":
    cli()
