import time

import pytest
from click.testing import CliRunner, Result

import aiida
import aiidalab_qe.__main__ as cli

# To learn more about testing click applications, see: https://click.palletsprojects.com/en/8.1.x/testing/


@pytest.mark.slow
def test_download_and_install_pseudos(tmp_path, aiida_profile, monkeypatch):
    """Test download pseudos to a specified directory and install them.
    And test install pseudos without download.
    Compare the time elapsed between the two methods.

    Act: Run the command `aiidalab-qe download-pseudos --dest <tmp_path>`.
    Assert: All pseudos are downloaded to the specified directory.

    Act: Run the command `aiidalab-qe install-pseudos --source <tmp_path>`.
    Assert: All pseudos are installed from the specified directory.

    Act: Run the command `aiidalab-qe install-pseudos`.
    Assert: All pseudos are installed from the default directory, but slow.

    Note: this test is slow, it takes about ~30 seconds to run.
    """
    from aiidalab_qe.setup.pseudos import (
        PSEUDODOJO_VERSION,
        SSSP_VERSION,
        pseudos_to_install,
    )

    # Mock the EXPECTED_PSEUDOS to speed up the test
    MOCK_EXPECTED_PSEUDOS = {
        f"SSSP/{SSSP_VERSION}/PBE/efficiency",
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/standard/upf",
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBE/SR/stringent/upf",
        f"PseudoDojo/{PSEUDODOJO_VERSION}/PBEsol/SR/stringent/upf",
    }
    # mock the EXPECTED_PSEUDOS
    monkeypatch.setattr(
        "aiidalab_qe.setup.pseudos.EXPECTED_PSEUDOS",
        MOCK_EXPECTED_PSEUDOS,
    )

    # clean the profile database
    aiida_profile.clear_profile()

    runner: CliRunner = CliRunner()

    start = time.time()
    result: Result = runner.invoke(cli.download_pseudos, ["--dest", str(tmp_path)])
    download_time = time.time() - start

    assert result.exit_code == 0
    assert "Pseudopotentials are downloaded!" in result.output

    files = [f for f in tmp_path.glob("**/*") if f.is_file()]

    for pseudo in MOCK_EXPECTED_PSEUDOS:
        filename = f"{pseudo.replace('/', '_')}.aiida_pseudo"
        assert filename in [f.name for f in files]

    # Install the pseudos from the tmp_path

    # Check that the pseudos are not installed yet
    assert len(pseudos_to_install()) == len(MOCK_EXPECTED_PSEUDOS)

    profile = aiida.get_profile()

    start = time.time()
    result: Result = runner.invoke(
        cli.install_pseudos, ["--source", tmp_path, "--profile", profile.name]
    )
    install_time = time.time() - start

    print(f"Download time: {download_time}, install time: {install_time}")

    assert result.exit_code == 0
    assert "Pseudopotentials are installed!" in result.output

    assert len(pseudos_to_install()) == 0

    # clean the profile database for the next test
    aiida_profile.clear_profile()

    # Check that the pseudos are not installed yet
    assert len(pseudos_to_install()) == len(MOCK_EXPECTED_PSEUDOS)

    profile = aiida.get_profile()
    start = time.time()
    result: Result = runner.invoke(cli.install_pseudos, ["--profile", profile.name])
    install_without_download_time = time.time() - start
    print(f"Install without download time: {install_without_download_time}")

    assert result.exit_code == 0
    assert "Pseudopotentials are installed!" in result.output

    assert len(pseudos_to_install()) == 0

    # Check the install without download is slower than install with download
    assert install_without_download_time > install_time
