import sys
import numpy as np
import unittest as ut
import espressomd
import espressomd.observables
from espressomd import utils
import tests_common

class TestCylindricalObservable(ut.TestCase):
    """
    Testcase for the cylindrical observables.

    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.time_step = 0.01
    system.box_l = [15.0, 15.0, 15.0]
    system.cell_system.skin = 0.4

    params = {
        'ids': range(100),
        'center': [7.5, 7.5, 7.5],  # center of the histogram
        'axis': 'y',
        'n_r_bins': 4,  # number of bins in r
        'n_phi_bins': 4,  # -*- in phi
        'n_z_bins': 4,  # -*- in z
        'min_r': 0.0,
        'min_phi': -np.pi,
        'min_z': -5.0,
        'max_r': 5.0,
        'max_phi': np.pi,
        'max_z': 5.0,
    }

    def swap_axis(self, arr, axis):
        if axis == 'x':
            arr = np.dot(tests_common.rotation_matrix([0, 1, 0], np.pi/2.0), arr)
        elif axis == 'y':
            arr = np.dot(tests_common.rotation_matrix([1, 0, 0], -np.pi/2.0), arr)
        return arr

    def swap_axis_inverse(self, arr, axis):
        if axis == 'x':
            arr = np.dot(tests_common.rotation_matrix([0, 1, 0], -np.pi/2.0), arr)
        elif axis == 'y':
            arr = np.dot(tests_common.rotation_matrix([1, 0, 0], np.pi/2.0), arr)
        return arr

    def pol_coords(self):
        positions = np.zeros((len(self.params['ids']), 3))
        for i, p in enumerate(self.system.part):
            tmp = p.pos - np.array(self.params['center'])
            tmp = self.swap_axis_inverse(tmp, self.params['axis'])
            positions[i, :] = tests_common.transform_pos_from_cartesian_to_polar_coordinates(tmp)
        return positions

    def set_particles(self):
        self.system.part.clear()
        # Parameters for an ellipse.
        a = 1.0  # semi minor-axis length
        b = 2.0  # semi major-axis length
        # Choose the cartesian velocities such that each particle gets the same
        # v_r, v_phi and v_z, respectively.
        self.v_r = .75
        self.v_phi = 2.5
        self.v_z = 1.5
        for i in range(len(self.params['ids'])):
            position = np.array([a *
                                 np.cos(i *
                                        2.0 *
                                        np.pi /
                                        (len(self.params['ids']) +
                                         1)), b *
                                 np.sin(i *
                                        2.0 *
                                        np.pi /
                                        (len(self.params['ids']) +
                                         1)), i *
                                  (self.params['max_z'] - self.params['min_z'])/
                                 (len(self.params['ids']) +
                                     1) - self.params['center'][2]])
            v_y = (position[0] * np.sqrt(position[0]**2.0 + position[1]**2.0) * \
                   self.v_phi + position[1] * self.v_r) / np.sqrt(position[0]**2.0 + position[1]**2.0)
            v_x = (self.v_r * np.sqrt(position[0]**2.0 + position[1] **
                                 2.0) - position[1] * v_y) / position[0]
            velocity = np.array([v_x, v_y, self.v_z])
            velocity = self.swap_axis(velocity, self.params['axis'])
            position = self.swap_axis(position, self.params['axis'])
            position += np.array(self.params['center'])
            self.system.part.add(id=i, pos=position, v=velocity)

    def calculate_numpy_histogram(self):
        pol_positions = self.pol_coords()
        np_hist, _ = np.histogramdd(pol_positions, bins=(self.params['n_r_bins'], self.params['n_phi_bins'], self.params['n_z_bins']),
                                    range=[(self.params['min_r'], self.params['max_r']), (self.params['min_phi'], self.params['max_phi']),
                                           (self.params['min_z'], self.params['max_z'])])
        return np_hist

    def normalize_with_bin_volume(self, histogram):
        bin_volume = tests_common.get_cylindrical_bin_volume(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'],
            self.params['min_r'],
            self.params['max_r'],
            self.params['min_phi'],
            self.params['max_phi'],
            self.params['min_z'],
            self.params['max_z'])
        for i in range(self.params['n_r_bins']):
            histogram[i, :, :] /= bin_volume[i]
        return histogram

    def density_profile_test(self):
        self.set_particles()
        # Set up the Observable.
        obs = espressomd.observables.CylindricalDensityProfile(**self.params)
        core_hist = np.array(
            obs.calculate()).reshape(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'])
        np_hist = self.calculate_numpy_histogram()
        np_hist = self.normalize_with_bin_volume(np_hist)
        np.testing.assert_array_almost_equal(np_hist, core_hist)
        self.assertEqual(obs.n_values(), len(np_hist.flatten()))

    def velocity_profile_test(self):
        self.set_particles()
        # Set up the Observable.
        obs = espressomd.observables.CylindricalVelocityProfile(**self.params)
        core_hist = np.array(
            obs.calculate()).reshape(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'],
            3)
        core_hist_v_r = core_hist[:, :, :, 0]
        core_hist_v_phi = core_hist[:, :, :, 1]
        core_hist_v_z = core_hist[:, :, :, 2]
        np_hist = self.calculate_numpy_histogram()
        for x in np.nditer(np_hist, op_flags=['readwrite']):
            if x[...] > 0.0:
                x[...] /= x[...]
        np.testing.assert_array_almost_equal(np_hist * self.v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(np_hist * self.v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(np_hist * self.v_z, core_hist_v_z)
        self.assertEqual(obs.n_values(), len(np_hist.flatten())*3)

    def flux_density_profile_test(self):
        self.set_particles()
        # Set up the Observable.
        obs = espressomd.observables.CylindricalFluxDensityProfile(**self.params)
        core_hist = np.array(
            obs.calculate()).reshape(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'],
            3)
        core_hist_v_r = core_hist[:, :, :, 0]
        core_hist_v_phi = core_hist[:, :, :, 1]
        core_hist_v_z = core_hist[:, :, :, 2]
        np_hist = self.calculate_numpy_histogram()
        np_hist = self.normalize_with_bin_volume(np_hist)
        np.testing.assert_array_almost_equal(np_hist * self.v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(np_hist * self.v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(np_hist * self.v_z, core_hist_v_z)
        self.assertEqual(obs.n_values(), len(np_hist.flatten())*3)

    def test_hist_x(self):
        self.params['axis'] = 'x'
        self.velocity_profile_test()
        self.flux_density_profile_test()
        self.density_profile_test()

    def test_hist_y(self):
        self.params['axis'] = 'y'
        self.velocity_profile_test()
        self.flux_density_profile_test()
        self.density_profile_test()

    def test_hist_z(self):
        self.params['axis'] = 'z'
        self.velocity_profile_test()
        self.flux_density_profile_test()
        self.density_profile_test()

if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        TestCylindricalObservable))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
