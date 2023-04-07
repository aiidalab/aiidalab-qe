import os
from pathlib import Path
from urllib.parse import urljoin

import pytest
import requests
import selenium.webdriver.support.expected_conditions as EC
from requests.exceptions import ConnectionError
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait


def is_responsive(url):
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return True
    except ConnectionError:
        return False


@pytest.fixture(scope="session")
def docker_compose_file(pytestconfig):
    return str(Path(pytestconfig.rootdir) / "tests_notebooks" / "docker-compose.yml")


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

    # install aiidalab-qe
    aiidalab_exec("pip install -U .", workdir=appdir)

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
        # By default, let's allow selenium functions to retry for 10s
        # till a given element is loaded, see:
        # https://selenium-python.readthedocs.io/waits.html#implicit-waits
        selenium.implicitly_wait(wait_time)
        window_width = 800
        window_height = 600
        selenium.set_window_size(window_width, window_height)

        selenium.find_element(By.ID, "ipython-main-app")
        selenium.find_element(By.ID, "notebook-container")
        WebDriverWait(selenium, 100).until(
            EC.invisibility_of_element((By.ID, "appmode-busy"))
        )

        return selenium

    return _selenium_driver


@pytest.fixture
def final_screenshot(request, screenshot_dir, selenium):
    """Take screenshot at the end of the test.
    Screenshot name is generated from the test function name
    by stripping the 'test_' prefix
    """
    screenshot_name = f"{request.function.__name__[5:]}.png"
    screenshot_path = Path.joinpath(screenshot_dir, screenshot_name)
    yield
    selenium.get_screenshot_as_file(screenshot_path)


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
