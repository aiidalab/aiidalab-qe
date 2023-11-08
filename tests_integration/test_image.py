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
    output = aiidalab_exec("verdi status", user=nb_user).decode().strip()
    assert "Connected to RabbitMQ" in output
    assert "Daemon is running" in output


@pytest.mark.usefixtures("notebook_service")
def test_pseudos_families_are_installed(aiidalab_exec, nb_user):
    # Check the aiida service is running and connected to RabbitMQ
    # The notebook_service fixture is needed to wait for the services to be up
    output = aiidalab_exec("aiida-pseudo list", user=nb_user).decode().strip()
    assert "SSSP" in output
    assert "PseudoDojo" in output

    # Two lines of header, 8 pseudos
    assert len(output.splitlines()) == 10
