"""For running the app from the command line used for post_install script.
"""

import click

from aiidalab_qe.app.common.setup_codes import codes_are_setup
from aiidalab_qe.app.common.setup_codes import install as install_qe_codes
from aiidalab_qe.app.submission.sssp import install as setup_sssp


@click.group()
def cli():
    pass


@cli.command()
@click.option("-f", "--force", is_flag=True)
def install_qe(force):
    try:
        for msg in install_qe_codes(force=force):
            click.echo(msg)
        assert codes_are_setup()
        click.secho("Codes are setup!", fg="green")
    except Exception as error:
        raise click.ClickException(f"Failed to set up QE failed: {error}")


@cli.command()
def install_sssp():
    try:
        for msg, _ in setup_sssp():
            click.echo(msg)
        click.secho("SSSP pseudo potentials are installed!", fg="green")
    except Exception as error:
        raise click.ClickException(f"Failed to set up pseudo potentials: {error}")


if __name__ == "__main__":
    cli()
