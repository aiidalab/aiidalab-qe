#!/usr/bin/env python
import time

from selenium.webdriver.common.by import By


def test_qe_app_take_screenshot(selenium, url):
    selenium.get(url("http://localhost:8100/apps/apps/quantum-espresso/qe.ipynb"))
    selenium.set_window_size(1920, 985)
    time.sleep(10)
    selenium.get_screenshot_as_file("screenshots/qe-app.png")


def test_qe_app_select_silicon(selenium, url):
    selenium.get(url("http://localhost:8100/apps/apps/quantum-espresso/qe.ipynb"))
    selenium.set_window_size(1920, 985)
    time.sleep(10)
    selenium.find_element(
        By.XPATH, '//li[@id="tab-key-17" and @class="lm-TabBar-tab p-TabBar-tab"]'
    ).click()  # click `From Examples` for input structure
    selenium.find_element(By.XPATH, "//option[@value='Diamond']").click()
    selenium.get_screenshot_as_file("screenshots/qe-app-select-diamond-selected.png")
    confirm_button = selenium.find_element(By.XPATH, "//button[contains(.,'Confirm')]")
    confirm_button.location_once_scrolled_into_view  # scroll into view
    confirm_button.click()
    selenium.get_screenshot_as_file("screenshots/qe-app-diamond-silicon-confirmed.png")
