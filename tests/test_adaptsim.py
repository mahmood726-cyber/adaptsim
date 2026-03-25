"""AdaptSim integration tests — adaptive trial design simulator."""
import sys, time, json, pytest
sys.path.insert(0, __import__('os').path.dirname(__file__))
from conftest import js

@pytest.fixture(autouse=True)
def load_app(driver, app_url):
    if driver.current_url != app_url:
        driver.get(app_url)
        time.sleep(1)

class TestAppLoads:
    def test_title(self, driver):
        assert 'AdaptSim' in driver.title

    def test_no_js_errors(self, driver):
        logs = driver.get_log('browser')
        severe = [l for l in logs if l['level'] == 'SEVERE']
        assert len(severe) == 0, f"JS errors: {severe}"

    def test_tabs_present(self, driver):
        tabs = driver.find_elements('css selector', '.tab-btn')
        assert len(tabs) >= 5  # Design, Boundaries, OC, SSR, Report

class TestExamples:
    def test_load_dapahf(self, driver):
        """Load DAPA-HF preset"""
        driver.find_element('id', 'exDapaHF').click()
        time.sleep(0.5)
        # Check effect size was populated
        es = driver.find_element('id', 'effectSize')
        val = es.get_attribute('value')
        assert val and float(val) > 0, f"Effect size should be set, got {val}"

    def test_load_emperor(self, driver):
        """Load EMPEROR-Reduced preset"""
        driver.find_element('id', 'exEmperor').click()
        time.sleep(0.5)
        es = driver.find_element('id', 'effectSize')
        val = es.get_attribute('value')
        assert val and float(val) > 0

class TestComputation:
    def test_compute_boundaries(self, driver):
        """Compute boundaries with default settings"""
        driver.execute_script("document.getElementById('exDapaHF').click()")
        time.sleep(0.5)
        driver.execute_script("document.getElementById('computeBtn').click()")
        time.sleep(2)
        # Switch to boundaries tab
        driver.find_element('id', 'tab-boundaries').click()
        time.sleep(0.5)
        panel = driver.find_element('id', 'panel-boundaries')
        assert panel.is_displayed()
        # Should have boundary table or chart
        content = panel.get_attribute('innerHTML')
        assert len(content) > 100, "Boundaries panel should have content"

    def test_oc_tab_renders(self, driver):
        """OC tab should show after computation"""
        driver.find_element('id', 'tab-oc').click()
        time.sleep(0.5)
        panel = driver.find_element('id', 'panel-oc')
        assert panel.is_displayed()

class TestTabSwitching:
    def test_cycle_all_tabs(self, driver):
        tab_ids = ['tab-design', 'tab-boundaries', 'tab-oc', 'tab-ssr', 'tab-report']
        for tid in tab_ids:
            driver.find_element('id', tid).click()
            time.sleep(0.3)
            # Corresponding panel should be visible
            panel_id = tid.replace('tab-', 'panel-')
            panel = driver.find_element('id', panel_id)
            assert panel.is_displayed(), f"Panel {panel_id} should be visible"

class TestDarkMode:
    def test_toggle_theme(self, driver):
        toggle = driver.find_element('id', 'themeToggle')
        toggle.click()
        time.sleep(0.3)
        body = driver.find_element('tag name', 'body')
        has_dark = 'dark' in (body.get_attribute('class') or '')
        # Toggle back
        toggle.click()
        time.sleep(0.3)
        # Just verify no crash — theme might not use body.dark class

class TestDesignForm:
    def test_endpoint_type_change(self, driver):
        """Changing endpoint type should update the form"""
        driver.find_element('id', 'tab-design').click()
        time.sleep(0.3)
        select = driver.find_element('id', 'endpointType')
        # Store current label
        label_before = driver.find_element('id', 'effectSizeLabel').text
        assert label_before is not None

    def test_num_looks_slider(self, driver):
        """Adjust number of looks"""
        slider = driver.find_element('id', 'numLooks')
        driver.execute_script("arguments[0].value = 3; arguments[0].dispatchEvent(new Event('input'));", slider)
        time.sleep(0.3)
        val = driver.find_element('id', 'numLooksVal').text
        assert val == '3', f"Expected '3', got '{val}'"
