import time
from pathlib import Path

import selenium.webdriver.support.expected_conditions as EC  # noqa: N812
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait


def test_qe_app_take_screenshot(selenium_driver, final_screenshot):
    driver = selenium_driver("qe.ipynb", wait_time=30.0)
    driver.set_window_size(1920, 1485)
    time.sleep(15)


def test_qe_app_select_silicon_and_confirm(
    selenium_driver,
    screenshot_dir,
    final_screenshot,
):
    driver = selenium_driver("qe.ipynb", wait_time=30.0)
    driver.set_window_size(1920, 1485)

    driver.find_element(By.CLASS_NAME, "qe-app-step-ready")  # ready on start

    # Open structure selection step
    element = WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.CLASS_NAME, "p-Accordion-child"))
    )
    element.click()
    # check that element has CSS class
    driver.find_element(By.CLASS_NAME, "qe-app-step-ready")  # still ready

    # Select the Silicon example
    element = WebDriverWait(driver, 60 * 2).until(
        EC.presence_of_element_located((By.XPATH, "//*[text()='From examples']"))
    )
    element.click()

    try:
        driver.find_element(By.XPATH, "//option[@value='Diamond']").click()
        time.sleep(10)
        # Selection configures the step
        driver.find_element(By.CLASS_NAME, "qe-app-step-configured")
        element = WebDriverWait(driver, 60).until(
            EC.element_to_be_clickable((By.XPATH, "//button[text()='Confirm']"))
        )
        assert element.is_enabled()
        element.click()
        # Confirming means step is successful
        driver.find_element(By.CLASS_NAME, "qe-app-step-success")
    except Exception:
        driver.find_element(By.TAG_NAME, "summary").click()
        time.sleep(10)
    # Scroll to bottom of screen
    driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
    # Take a screenshot of the selected diamond
    driver.get_screenshot_as_file(
        str(Path.joinpath(screenshot_dir, "qe-app-select-diamond-selected.png"))
    )
    # Test that we have indeed proceeded to the next step
    driver.find_element(By.XPATH, "//span[contains(.,'✓ Step 1')]")
