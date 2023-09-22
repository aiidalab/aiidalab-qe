import requests


def test_notebook_service_available(notebook_service):
    url, token = notebook_service
    response = requests.get(f"{url}/?token={token}")
    assert response.status_code == 200


def test_verdi_status(aiidalab_exec, nb_user, notebook_service):
    # Check the aiida service is running and connected to RabbitMQ
    # The notebook_service fixture is needed to wait for the services to be up
    output = aiidalab_exec("verdi status", user=nb_user).decode().strip()
    assert "Connected to RabbitMQ" in output
    assert "Daemon is running" in output
