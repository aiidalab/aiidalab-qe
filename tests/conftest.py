import os
from pathlib import Path
from urllib.parse import urljoin

import pytest
import requests
from requests.exceptions import ConnectionError
from selenium.webdriver.common.by import By


def is_responsive(url):
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return True
    except ConnectionError:
        return False


@pytest.fixture(scope="session")
def docker_compose(docker_services):
    return docker_services._docker_compose


@pytest.fixture(scope="session")
def aiidalab_exec(docker_compose):
    def execute(command, user=None, workdir=None, **kwargs):
        opts = "-T"
        if user:
            opts = f"{opts} --user={user}"
        if workdir:
            opts = f"{opts} --workdir={workdir}"
        command = f"exec {opts} aiidalab {command}"

        return docker_compose.execute(command, **kwargs)

    return execute


@pytest.fixture(scope="session")
def nb_user(aiidalab_exec):
    return aiidalab_exec("bash -c 'echo \"${NB_USER}\"'").decode().strip()


@pytest.fixture(scope="session")
def notebook_service(docker_ip, docker_services, aiidalab_exec, nb_user):
    """Ensure that HTTP service is up and responsive."""

    # Directory ~/apps/aiidalab-qe/ is mounted by docker,
    # make it writeable for jovyan user, needed for `pip install`
    appdir = f"/home/{nb_user}/apps/aiidalab-qe"
    aiidalab_exec(f"chmod -R a+rw {appdir}", user="root")

    # Install workchains
    aiidalab_exec("pip install .", workdir=f"{appdir}/src", user=nb_user)

    # Install App
    install_command = "bash -c 'python tests/helper_dep_requirements.py && pip install -r /tmp/requirements.txt'"
    aiidalab_exec(install_command, workdir=appdir, user=nb_user)

    # `port_for` takes a container port and returns the corresponding host port
    port = docker_services.port_for("aiidalab", 8888)
    url = f"http://{docker_ip}:{port}"
    token = os.environ["JUPYTER_TOKEN"]
    docker_services.wait_until_responsive(
        timeout=30.0, pause=0.1, check=lambda: is_responsive(url)
    )
    return url, token


@pytest.fixture(scope="function")
def selenium_driver(selenium, notebook_service):
    def _selenium_driver(nb_path, wait_time=5.0):
        url, token = notebook_service
        url_with_token = urljoin(url, f"apps/apps/aiidalab-qe/{nb_path}?token={token}")
        selenium.get(f"{url_with_token}")
        selenium.implicitly_wait(wait_time)  # must wait until the app loaded

        selenium.find_element(By.ID, "ipython-main-app")
        selenium.find_element(By.ID, "notebook-container")

        return selenium

    return _selenium_driver


@pytest.fixture(scope="session")
def screenshot_dir():
    sdir = Path.joinpath(Path.cwd(), "screenshots")
    try:
        os.mkdir(sdir)
    except FileExistsError:
        pass
    return sdir


@pytest.fixture
def firefox_options(firefox_options):
    firefox_options.add_argument("--headless")
    return firefox_options


@pytest.fixture
def chrome_options(chrome_options):
    chrome_options.add_argument("--headless")
    return chrome_options
