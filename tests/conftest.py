"""Shared Selenium fixtures for AdaptSim tests."""
import os

import pytest
from selenium import webdriver
from selenium.common.exceptions import WebDriverException
from selenium.webdriver.chrome.options import Options

HTML = os.path.join(os.path.dirname(__file__), '..', 'adaptsim.html')

@pytest.fixture(scope='session')
def driver():
    opts = Options()
    opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--disable-gpu')
    opts.set_capability('goog:loggingPrefs', {'browser': 'ALL'})
    try:
        d = webdriver.Chrome(options=opts)
    except WebDriverException as exc:
        pytest.skip(f"Chrome WebDriver unavailable: {exc}")
    d.implicitly_wait(3)
    yield d
    d.quit()

@pytest.fixture(scope='session')
def app_url():
    return 'file:///' + os.path.abspath(HTML).replace('\\', '/')


def js(driver, script):
    return driver.execute_script('return ' + script)
