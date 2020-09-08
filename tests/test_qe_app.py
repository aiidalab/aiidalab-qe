#!/usr/bin/env python
from selenium.webdriver.common.by import By
import time


def test_qe_app_select_silicon(selenium, url):
    selenium.get(url("http://localhost:8100/apps/apps/quantum-espresso/qe.ipynb"))
    selenium.set_window_size(1920, 985)
    selenium.find_element(By.CSS_SELECTOR, ".p-TabBar-tab:nth-child(6) > .p-TabBar-tabLabel").click()
    selenium.find_element(By.XPATH, "//option[@value=\'Silicon\']").click()
    selenium.find_element(By.XPATH, "//div[@id=\'notebook-container\']/div/div[2]/div[2]/div[2]/div[3]/div/div[2]/div/div[2]/div/button").click()
    selenium.find_element(By.XPATH, "//div[@id=\'notebook-container\']/div/div[2]/div[2]/div[2]/div[3]/div/div[2]/div[2]/div[2]/div/div[2]/div/div[2]/div/div/ul/li[3]/div[2]").click()
    selenium.get_screenshot_as_file('screenshots/qe-app-select-silicon.png')
