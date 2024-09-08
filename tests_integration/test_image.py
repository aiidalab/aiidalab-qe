import pytest
import requests


def test_notebook_service_available(notebook_service):
    url, token = notebook_service
    response = requests.get(f"{url}/?token={token}")
    assert response.status_code == 200


@pytest.mark.usefixtures("notebook_service")
def test_verdi_status(aiidalab_exec, nb_user):
    # Check the aiida service is running and connected to RabbitMQ
    # The notebook_service fixture is needed to wait for the services to be up
    output = aiidalab_exec("verdi status", user=nb_user)
    for status in ("version", "config", "profile", "storage", "broker", "daemon"):
        assert f"✔ {status}" in output
    assert "/home/jovyan/.aiida" in output
    assert "Daemon is running" in output
    assert "Unable to connect to broker" not in output


@pytest.mark.usefixtures("notebook_service")
def test_pseudos_families_are_installed(aiidalab_exec, nb_user):
    # Check the aiida service is running and connected to RabbitMQ
    # The notebook_service fixture is needed to wait for the services to be up
    output = aiidalab_exec("aiida-pseudo list", user=nb_user)
    assert "SSSP" in output
    assert "PseudoDojo" in output

    # Two lines of header, 8 pseudos
    assert len(output.splitlines()) == 14
