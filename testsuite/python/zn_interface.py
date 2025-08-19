import unittest
import time
import espressomd.zn
from selenium import webdriver
from selenium.webdriver.common.by import By

system = espressomd.System(box_l=[1]*3)
vis = espressomd.zn.Visualizer(system)
url = vis.address

#Default Methods
method_list = [['ConnectedParticles', 'NoneSelection', 'All', 'Invert', 'Range', 'Random', 'IdenticalSpecies', 'Neighbour'],
               ['Delete', 'Rotate', 'Translate', 'Duplicate', 'ChangeType', 'AddLineParticles',
                   'Wrap', 'Center', 'Replicate', 'Connect', 'NewCanvas', 'RemoveAtoms'],
               ['Plane', 'Sphere', 'Box', 'Circle', 'Cone', 'Cylinder', 'Dodecahedron', 'Icosahedron',
                   'Octahedron', 'Ring', 'Tetrahedron', 'Torus', 'TorusKnot', 'Rhomboid', 'Ellipsoid'],
               ['Properties1D', 'DihedralAngle', 'Distance', 'Properties2D']]


def s(t=1):
    time.sleep(t)


class RegisterTestCase(unittest.TestCase):
    def setUp(self):
        chrome_options = webdriver.ChromeOptions()
        chrome_options.add_argument("headless") 
        self.driver = webdriver.Chrome(options=chrome_options)
        self.driver.get(url)

    def tearDown(self):
        self.driver.quit()

    def test_register(self):
        next_option_names = []
        idx = 0
        driver = self.driver
        buttons = driver.find_elements(
            By.XPATH, '//button[@class="btn btn-outline-tertiary"]')
        for button in buttons:
            button.click()
            s(t=0.1)
            elements = driver.find_elements(By.ID, "root[method]switcher")
            j = len(elements) - 1
            options = [x for x in elements[j].find_elements(
                By.TAG_NAME, "option")]
            option_names = [x.get_attribute("value") for x in options]
            if option_names == next_option_names:
                continue
            else:
                assert option_names == method_list[idx], f"Difference: {
                    set(option_names) ^ set(method_list[idx])}"
                idx += 1
                next_option_names = option_names


if __name__ == "__main__":
    unittest.main()
