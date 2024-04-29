import time
from pathlib import Path

import selenium.webdriver.support.expected_conditions as EC
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

    element = WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.XPATH, "//*[text()='From Examples']"))
    )
    element.click()

    try:
        driver.find_element(
            By.XPATH, "//option[@value='Diamond (primitive cell)']"
        ).click()
        time.sleep(10)
        element = WebDriverWait(driver, 60).until(
            EC.element_to_be_clickable((By.XPATH, "//button[text()='Confirm']"))
        )
        element.click()
    except Exception:
        driver.find_element(By.TAG_NAME, "summary").click()
        time.sleep(10)
    # Take a screenshot of the selected diamond
    driver.get_screenshot_as_file(
        str(Path.joinpath(screenshot_dir, "qe-app-select-diamond-selected.png"))
    )
    # Test that we have indeed proceeded to the next step
    driver.find_element(By.XPATH, "//span[contains(.,'✓ Step 1')]")
