"""For running the app from the command line used for post_install script.
"""

import click
from aiida import load_profile

from aiidalab_qe.common.setup_codes import codes_are_setup
from aiidalab_qe.common.setup_codes import install as install_qe_codes
from aiidalab_qe.common.setup_pseudos import install as setup_pseudos


@click.group()
def cli():
    pass


@cli.command()
@click.option("-f", "--force", is_flag=True)
@click.option("-p", "--profile", default="default")
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
@click.option("-p", "--profile", default="default")
def install_pseudos(profile):
    load_profile(profile)
    try:
        for msg, _ in setup_pseudos():
            click.echo(msg)
        click.secho("Pseudopotentials are installed!", fg="green")
    except Exception as error:
        raise click.ClickException(f"Failed to set up pseudo potentials: {error}")


if __name__ == "__main__":
    cli()
