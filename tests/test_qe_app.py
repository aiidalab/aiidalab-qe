from pathlib import Path
import requests
from selenium.webdriver.common.by import By


def test_notebook_service_available(notebook_service):
    url, token = notebook_service
    response = requests.get(f"{url}/?token={token}")
    assert response.status_code == 200


def test_qe_app_take_screenshot(selenium_driver, screenshot_dir):
    driver = selenium_driver("qe.ipynb", wait_time=30.0)
    driver.set_window_size(1920, 985)
    driver.get_screenshot_as_file(Path.joinpath(screenshot_dir, "qe-app.png"))


def test_qe_app_select_silicon(selenium_driver, screenshot_dir):
    driver = selenium_driver("qe.ipynb", wait_time=30.0)
    driver.set_window_size(1920, 985)
    driver.find_element(
        By.XPATH, "//*[text()='From Examples']"
    ).click()  # click `From Examples` tab for input structure
    driver.find_element(By.XPATH, "//option[@value='Diamond']").click()
    driver.get_screenshot_as_file(Path.joinpath(screenshot_dir, "qe-app-select-diamond-selected.png"))
    confirm_button = driver.find_element(By.XPATH, "//button[text()='Confirm']")
    confirm_button.location_once_scrolled_into_view  # scroll into view
    confirm_button.click()
    driver.get_screenshot_as_file(Path.joinpath(screenshot_dir, "qe-app-select-diamond-confirmed.png"))
