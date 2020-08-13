#!/usr/bin/env python
from selenium.webdriver.common.by import By
import time


def test_qe_app_select_silicon(selenium, url):
    selenium.get(url("http://localhost:8100/apps/apps/quantum-espresso/qe.ipynb"))
    selenium.set_window_size(1200, 941)
    selenium.find_element(By.CSS_SELECTOR, ".p-TabBar-tab:nth-child(5) > .p-TabBar-tabLabel").click()
    selenium.find_element(By.XPATH, "//div[@id=\'notebook-container\']/div/div[2]/div[2]/div[2]/div[3]/div/div[2]/div/div[2]/div/div[2]/div/div[2]/div[5]/div/select").click()
    dropdown = selenium.find_element(By.XPATH, "//div[@id=\'notebook-container\']/div/div[2]/div[2]/div[2]/div[3]/div/div[2]/div/div[2]/div/div[2]/div/div[2]/div[5]/div/select")
    dropdown.find_element(By.XPATH, "//option[. = 'Silicon']").click()
    time.sleep(1)
    selenium.get_screenshot_as_file('screenshots/qe-app-select-silicon.png')
    selenium.find_element(By.CSS_SELECTOR, ".widget-button:nth-child(4)").click()
