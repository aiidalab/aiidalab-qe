import pytest
from click.testing import CliRunner

import aiida
import aiidalab_qe.__main__ as cli
from aiida.manage.configuration import get_profile

# To learn more about testing click applications, see: https://click.palletsprojects.com/en/8.1.x/testing/


def _mock_expected_pseudos(monkeypatch):
    from aiidalab_qe.setup.pseudos import (
        PSEUDODOJO_VERSION,
        SSSP_VERSION,
    )

    expected_pseudos = {
        f"SSSP/{SSSP_VERSION}/PBE/efficiency",
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf",
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/stringent/upf",
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/stringent/upf",
    }
    monkeypatch.setattr(
        "aiidalab_qe.setup.pseudos.EXPECTED_PSEUDOS",
        expected_pseudos,
    )
    return expected_pseudos


@pytest.mark.slow
@pytest.mark.usefixtures("aiida_profile_clean")
def test_download_and_install_pseudos_from_directory(tmp_path, monkeypatch):
    """Download pseudo archives and install them from the temporary directory."""
    from aiidalab_qe.setup.pseudos import pseudos_to_install

    expected_pseudos = _mock_expected_pseudos(monkeypatch)
    runner = CliRunner()

    result = runner.invoke(cli.download_pseudos, ["--dest", str(tmp_path)])

    assert result.exit_code == 0
    assert "Pseudopotentials are downloaded!" in result.output

    files = [f for f in tmp_path.glob("**/*") if f.is_file()]

    for pseudo in expected_pseudos:
        filename = f"{pseudo.replace('/', '_')}.aiida_pseudo"
        assert filename in [f.name for f in files]

    # Install the pseudos from the tmp_path

    # Check that the pseudos are not installed yet
    assert len(pseudos_to_install()) == len(expected_pseudos)

    profile = get_profile()
    result = runner.invoke(
        cli.install_pseudos, ["--source", tmp_path, "--profile", profile.name]
    )

    assert result.exit_code == 0
    assert "Pseudopotentials are installed!" in result.output

    assert len(pseudos_to_install()) == 0


@pytest.mark.slow
@pytest.mark.usefixtures("aiida_profile_clean")
def test_install_pseudos_from_default_source(monkeypatch):
    """Install pseudos from the packaged default source into a clean profile."""
    from aiidalab_qe.setup.pseudos import pseudos_to_install

    expected_pseudos = _mock_expected_pseudos(monkeypatch)
    assert len(pseudos_to_install()) == len(expected_pseudos)
    profile = get_profile()
    result = CliRunner().invoke(cli.install_pseudos, ["--profile", profile.name])

    assert result.exit_code == 0
    assert "Pseudopotentials are installed!" in result.output

    assert len(pseudos_to_install()) == 0


@pytest.mark.parametrize("post_install", ["setup", None])
def test_test_plugin_runs_post_install_verifier_when_configured(
    monkeypatch, post_install
):
    from types import SimpleNamespace

    from aiidalab_qe.plugins import registry
    from aiidalab_qe.plugins import utils as plugin_utils

    calls = []
    monkeypatch.setattr(aiida, "load_profile", lambda _profile: None)
    monkeypatch.setattr(
        plugin_utils, "test_plugin_functionality", lambda _name: (True, "")
    )
    monkeypatch.setattr(
        registry,
        "get_plugin_registry",
        lambda: SimpleNamespace(
            plugins={"my-plugin": SimpleNamespace(post_install=post_install)}
        ),
    )
    monkeypatch.setattr(
        cli.subprocess,
        "run",
        lambda command, check: calls.append((command, check)),
    )

    result = CliRunner().invoke(cli.test_plugin, ["my-plugin", "--profile", "test"])

    assert result.exit_code == 0
    if post_install:
        assert calls == [
            (
                [
                    cli.sys.executable,
                    "-m",
                    "my_plugin",
                    "verify-post-install",
                    "--profile",
                    "test",
                ],
                True,
            )
        ]
    else:
        assert calls == []
