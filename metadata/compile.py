#!/usr/bin/env python
import json
import re
import shutil
from pathlib import Path
from subprocess import run, CalledProcessError

import click
from dulwich.repo import Repo
from dulwich.porcelain import active_branch, reset
from jinja2 import Environment, FileSystemLoader
from jsonschema import validate
from packaging.specifiers import SpecifierSet
from packaging.version import parse

env = Environment()



def get_version_identifier(version):
    "Get version identifier from version (e.g. refs/tags/v1.0.0 -> v1.0.0)."
    if version.startswith('refs/tags/'):
        return version[len('refs/tags/'):]
    if version.startswith('refs/heads/'):
        return version[len('refs/heads/'):]
    if version.startswith('refs/remotes/'):  # remote branch
        return re.sub(r'refs\/remotes\/(.+?)\/', '', version)
    return version


def get_requirements_for_ref(repo, ref):
    """Get the requirements (from requirements.txt) for a given ref."""
    try:
        return tuple(
            run(['git', 'show', f'{ref.decode()}:requirements.txt'],
                check=True, capture_output=True, encoding='utf-8').stdout.splitlines())
    except CalledProcessError as error:
        if error.returncode == 128:
            return tuple()  # file does not exist for that ref
        raise


@click.command()
@click.option('-m', '--metadata-template', default='metadata.json.in')
def cli(metadata_template):
    """Generate the 'metadata.json' file for this app."""
    metadata_dir = Path(__file__).resolve().parent
    output_dir = metadata_dir / 'build'
    root = metadata_dir.parent

    repo = Repo(root)
    versions = {parse(get_version_identifier(ref.decode())): ref for ref in repo.get_refs()}
    requirements = {version: get_requirements_for_ref(repo, ref) for version, ref in versions.items()}

    def get_requirements(spec, version=None):
        spec = SpecifierSet(spec)
        if version is None:
            matching_versions = [version for version in sorted(versions) if version in spec]
            matching_requirements = {requirements[version] for version in matching_versions}
            if len(matching_requirements) == 0:
                raise RuntimeError(f"Unable to determine requirements for specifier '{spec}'.")
            elif len(matching_requirements) > 1:
                raise RuntimeError(f"Requirements for specifier '{spec}' are not uniform.")
            reqs = matching_requirements.pop()
        else:
            reqs = requirements[parse(version)]

        return json.dumps({str(spec): reqs})[1:-1]


    env = Environment(loader=FileSystemLoader(metadata_dir / 'templates'))
    env.filters['get_requirements'] = get_requirements

    # Writing output...
    output_dir.mkdir(exist_ok=True)

    # index.html
    index_html = output_dir / 'index.html'
    index_html.write_text(env.get_template('index.html').render())
    click.echo(f"Write {index_html.relative_to(Path.cwd())}")

    # metadata.json
    metadata_json = output_dir / 'metadata.json'
    metadata_template = env.get_template(metadata_template)
    try:
        metadata = json.loads(metadata_template.render())
    except json.decoder.JSONDecodeError as error:
        raise RuntimeError(f"{error}\n{metadata_template.render()}")
    validate(
        instance=metadata,
        schema={"$ref": "https://aiidalab.github.io/aiidalab-registry/schemas/v1/metadata.schema.json"})
    metadata_json.write_text(json.dumps(metadata, indent=2))
    click.echo(f"Write {metadata_json.relative_to(Path.cwd())}")

    # QE.jpg
    src = root / 'miscellaneous' / 'logos' / 'QE.jpg'
    dst = output_dir / 'miscellaneous' / 'logos' / 'QE.jpg'
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    click.echo(f"Copy {dst.relative_to(Path.cwd())}")


if __name__ == '__main__':
    cli()
