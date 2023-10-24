import aiida
import pytest
from click.testing import CliRunner, Result

import aiidalab_qe.__main__ as cli

# To learn more about testing click applications, see: https://click.palletsprojects.com/en/8.1.x/testing/


@pytest.mark.usefixtures("aiida_profile_clean")
def test_download_and_install_pseudos(tmp_path):
    """Test download pseudos to a specified directory.

    Act: Run the command `aiidalab-qe download-pseudos --dest <tmp_path>`.
    Assert: All pseudos are downloaded to the specified directory.

    Act: Run the command `aiidalab-qe install-pseudos --source <tmp_path>`.
    Assert: All pseudos are installed from the specified directory.
    """
    from aiidalab_qe.common.setup_pseudos import EXPECTED_PSEUDOS, pseudos_to_install

    runner: CliRunner = CliRunner()
    result: Result = runner.invoke(cli.download_pseudos, ["--dest", str(tmp_path)])
    assert result.exit_code == 0
    assert "Pseudopotentials are downloaded!" in result.output

    files = [f for f in tmp_path.glob("**/*") if f.is_file()]

    for pseudo in EXPECTED_PSEUDOS:
        filename = f"{pseudo.replace('/', '_')}.aiida_pseudo"
        assert filename in [f.name for f in files]

    # Install the pseudos from the tmp_path

    # Check that the pseudos are not installed yet
    assert len(pseudos_to_install()) == 8

    profile = aiida.get_profile()
    result: Result = runner.invoke(
        cli.install_pseudos, ["--source", tmp_path, "--profile", profile.name]
    )
    assert result.exit_code == 0
    assert "Pseudopotentials are installed!" in result.output

    assert len(pseudos_to_install()) == 0
