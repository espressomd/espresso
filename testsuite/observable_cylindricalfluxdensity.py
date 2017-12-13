import sys
import numpy as np
import unittest as ut
import espressomd
import espressomd.observables
from espressomd import utils
import tests_common

class TestCylindricalFluxDensityObservable(ut.TestCase):
    """
    Testcase for the CylindricalFluxDensityObservable.

    """
    system = espressomd.System()
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    system.time_step = 0.01
    system.box_l = [15.0, 15.0, 15.0]
    system.cell_system.skin = 0.4

    params = {
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
        'N': 100  # number of particles
    }


    def swap_axis(self, arr, axis):
        if axis == 'x':
            arr = np.array([arr[2], arr[1], -arr[0]])
        elif axis == 'y':
            arr = np.array([arr[0], arr[2], -arr[1]])
        return arr

    def swap_axis_inverse(self, arr, axis):
        if axis == 'x':
            arr = np.array([-arr[2], arr[1], arr[0]])
        elif axis == 'y':
            arr = np.array([arr[0], -arr[2], arr[1]])
        return arr


    def pol_coords(self):
        positions = np.zeros((self.params['N'], 3))
        velocities = np.zeros((self.params['N'], 3))
        for i, p in enumerate(self.system.part):
            tmp = p.pos - np.array(self.params['center'])
            tmp = self.swap_axis_inverse(tmp, self.params['axis'])
            positions[i, :] = tests_common.transform_pos_from_cartesian_to_polar_coordinates(tmp)
        return positions


    def hist_test(self):
        # Calculate the histogram normalization.
        bin_volume = np.zeros(self.params['n_r_bins'])
        r_bin_size = (
            self.params['max_r'] - self.params['min_r']) / self.params['n_r_bins']
        phi_bin_size = (
            self.params['max_phi'] - self.params['min_phi']) / self.params['n_phi_bins']
        z_bin_size = (
            self.params['max_z'] - self.params['min_z']) / self.params['n_z_bins']
        for i in range(self.params['n_r_bins']):
            bin_volume[i] = np.pi * ((self.params['min_r'] + r_bin_size * (i + 1))**2.0 -
                                     (self.params['min_r'] + r_bin_size * i)**2.0) * \
                phi_bin_size / (2.0 * np.pi) * z_bin_size
        self.system.part.clear()
        # Parameters for an ellipse.
        a = 1.0  # semi minor-axis length
        b = 2.0  # semi major-axis length
        # Choose the cartesian velocities such that each particle gets the same
        # v_r, v_phi and v_z, respectively.
        v_r = .75
        v_phi = 2.5
        v_z = 1.5
        for i in range(self.params['N']):
            position = np.array([a *
                                 np.cos(i *
                                        2.0 *
                                        np.pi /
                                        (self.params['N'] +
                                         1)), b *
                                 np.sin(i *
                                        2.0 *
                                        np.pi /
                                        (self.params['N'] +
                                         1)), i *
                                  (self.params['max_z'] - self.params['min_z'])/
                                 (self.params['N'] +
                                     1) - self.params['center'][2]])
            v_y = (position[0] * np.sqrt(position[0]**2.0 + position[1]**2.0) * \
                   v_phi + position[1] * v_r) / np.sqrt(position[0]**2.0 + position[1]**2.0)
            v_x = (v_r * np.sqrt(position[0]**2.0 + position[1] **
                                 2.0) - position[1] * v_y) / position[0]
            velocity = np.array([v_x, v_y, v_z])
            velocity = self.swap_axis(velocity, self.params['axis'])
            position = self.swap_axis(position, self.params['axis'])
            position += np.array(self.params['center'])
            self.system.part.add(id=i, pos=position, v=velocity)
        pol_positions = self.pol_coords()
        # Set up the Observable.
        p = espressomd.observables.CylindricalFluxDensityProfile(
            ids=range(
                self.params['N']),
            center=self.params['center'],
            axis=self.params['axis'],
            n_r_bins=self.params['n_r_bins'],
            n_phi_bins=self.params['n_phi_bins'],
            n_z_bins=self.params['n_z_bins'],
            min_r=self.params['min_r'],
            min_phi=self.params['min_phi'],
            min_z=self.params['min_z'],
            max_r=self.params['max_r'],
            max_phi=self.params['max_phi'],
            max_z=self.params['max_z'])
        core_hist = np.array(
            p.calculate()).reshape(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'],
            3)
        core_hist_v_r = core_hist[:, :, :, 0]
        core_hist_v_phi = core_hist[:, :, :, 1]
        core_hist_v_z = core_hist[:, :, :, 2]

        np_hist, _ = np.histogramdd(pol_positions, bins=(self.params['n_r_bins'], self.params['n_phi_bins'], self.params['n_z_bins']),
                                    range=[(self.params['min_r'], self.params['max_r']), (self.params['min_phi'], self.params['max_phi']),
                                           (self.params['min_z'], self.params['max_z'])])
        # Normalization
        for i in range(self.params['n_r_bins']):
            np_hist[i, :, :] /= bin_volume[i]
        np.testing.assert_array_almost_equal(np_hist * v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(np_hist * v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(np_hist * v_z, core_hist_v_z)

    def test_hist_x(self):
        self.params['axis'] = 'x'
        self.hist_test()

    def test_hist_y(self):
        self.params['axis'] = 'y'
        self.hist_test()

    def test_hist_z(self):
        self.params['axis'] = 'z'
        self.hist_test()

if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        TestCylindricalFluxDensityObservable))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
