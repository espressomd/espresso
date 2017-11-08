import sys
import numpy as np
import unittest as ut
import espressomd
import espressomd.observables
from espressomd import utils


class TestCylindricalFluxDensityObservable(ut.TestCase):
    """
    Testcase for the CylindricalFluxDensityObservable.

    """
    system = espressomd.System()
    system.time_step = 0.01
    system.box_l = [10.0, 10.0, 10.0]
    system.cell_system.skin = 0.4

    params = {
        'center': [5.0, 5.0, 0.0],  # center of the histogram
        'n_r_bins': 10,  # number of bins in r
        'n_phi_bins': 1,  # -*- in phi
        'n_z_bins': 1,  # -*- in z
        'min_r': 0.0,
        'min_phi': -np.pi,
        'min_z': 0.0,
        'max_r': 5.0,
        'max_phi': np.pi,
        'max_z': 10.0,
        'N': 100  # number of particles
    }

    def get_ind(self, n_r_bins, n_phi_bins, n_z_bins, ind):
        return utils.get_unravelled_index(
            [n_r_bins, n_phi_bins, n_z_bins, 3], 4, ind)

    def test_r(self):
        """Test binning in r.

        """
        self.system.part.clear()
        for i in range(self.params['N']):
            self.system.part.add(id=i, pos=np.array(
                self.params['center']) + np.array([i * 5.0 / self.params['N'], 0, 0]), v=[0.0, 0.0, 3.0])
        # Set up the Observable.
        self.p = espressomd.observables.CylindricalFluxDensityProfile(
            ids=range(
                self.params['N']),
            center=self.params['center'],
            n_r_bins=self.params['n_r_bins'],
            n_phi_bins=self.params['n_phi_bins'],
            n_z_bins=self.params['n_z_bins'],
            min_r=self.params['min_r'],
            min_phi=self.params['min_phi'],
            min_z=self.params['min_z'],
            max_r=self.params['max_r'],
            max_phi=self.params['max_phi'],
            max_z=self.params['max_z'])

        self.r_bin_size = (
            self.params['max_r'] - self.params['min_r']) / self.params['n_r_bins']
        self.hist = np.array(self.p.calculate())
        self.res = []
        for ind, v_dens in enumerate(self.hist):
            # Get the unflattened index.
            _ind = self.get_ind(self.params['n_r_bins'], self.params['n_phi_bins'], self.params['n_z_bins'], ind)
            # If we have the 'z' component add the velocity density to the
            # result.
            if _ind[3] == 2:
                # Append the mid of the bin position and the z-component of the
                # velocity density.
                self.res.append([(_ind[0] + 0.5) * self.r_bin_size, v_dens])
        self.res = np.array(self.res)
        # The Riemann sum should yield the sum of all z-velocities
        riemann_sum = np.sum(self.res[:, 0] * self.res[:, 1]) * 2.0 * np.pi * \
            (self.params['max_z'] - self.params['min_z']) * self.r_bin_size
        self.assertAlmostEqual(riemann_sum, self.params['N'] * 3.0, places=4)
        # Despite the Jacobian (r) factor all values should be the same by construction.
        tmp = self.res[:, 0] * self.res[:, 1]
        self.assertTrue(all(x - tmp[0] < 1E-8 for x in tmp))

    def test_phi(self):
        """Test binning in phi.

        """
        self.system.part.clear()
        # Place particles on a circle around the center.
        angles = np.linspace(0.0, 2.0*np.pi, num=self.params['N'], endpoint=False)
        for i in range(self.params['N']):
            pos = np.array(self.params['center'])
            pos += np.array([np.cos(angles[i]),
                             np.sin(angles[i]), 0.0])
            self.system.part.add(id=i, pos=pos, v=[0.0, 0.0, 3.0])
        # Set up the Observable.
        p = espressomd.observables.CylindricalFluxDensityProfile(
            ids=range(
                self.params['N']),
            center=self.params['center'],
            n_r_bins=1,
            n_phi_bins=10,
            n_z_bins=self.params['n_z_bins'],
            min_r=self.params['min_r'],
            min_phi=self.params['min_phi'],
            min_z=self.params['min_z'],
            max_r=self.params['max_r'],
            max_phi=self.params['max_phi'],
            max_z=self.params['max_z'])
        self.phi_bin_size = (
            self.params['max_phi'] - self.params['min_phi']) / 10.0
        self.hist = np.array(p.calculate())
        self.res = []
        for ind, v_dens in enumerate(self.hist):
            _ind = self.get_ind(1, 10, self.params['n_z_bins'], ind)
            if _ind[3] == 2:
                print _ind, v_dens
                self.res.append([(_ind[1] + 0.5) * self.phi_bin_size, v_dens])
        self.res = np.array(self.res)
        # The Riemann sum should yield the sum of all z-velocities.
        riemann_sum = np.sum(self.res[:, 1]) * \
            0.5 * 5.0 * 5.0 * 10.0 * self.phi_bin_size
        self.assertAlmostEqual(riemann_sum, self.params['N'] * 3.0, places=4)
        # All values should be the same by construction.
        #self.assertTrue(all(x - self.res[0, 1] < 1E-8 for x in self.res[:, 1]))

    def test_z(self):
        """Test binning in z.

        """
        self.system.part.clear()
        for i in range(self.params['N']):
            pos = np.array(self.params['center'])
            pos += np.array([0.0, 0.0, i * self.params['max_z'] / self.params['N']])
            self.system.part.add(id=i, pos=pos, v=[0.0, 0.0, 3.0])
        # Set up the Observable.
        p = espressomd.observables.CylindricalFluxDensityProfile(
            ids=range(
                self.params['N']),
            center=self.params['center'],
            n_r_bins=1,
            n_phi_bins=1,
            n_z_bins=10,
            min_r=self.params['min_r'],
            min_phi=self.params['min_phi'],
            min_z=self.params['min_z'],
            max_r=self.params['max_r'],
            max_phi=self.params['max_phi'],
            max_z=self.params['max_z'])
        self.z_bin_size = (
            self.params['max_z'] - self.params['min_z']) / 10.0
        self.hist = np.array(p.calculate())
        self.res = []
        for ind, v_dens in enumerate(self.hist):
            _ind = self.get_ind(1, 1, 10, ind)
            if _ind[3] == 2:
                self.res.append([(_ind[2] + 0.5) * self.z_bin_size, v_dens])
        self.res = np.array(self.res)
        # The Riemann sum should yield the sum of all z-velocities
        riemann_sum = np.sum(self.res[:, 1]) * 2.0 * \
            np.pi * 0.5 * 5.0 * 5.0 * self.z_bin_size
        self.assertAlmostEqual(riemann_sum, self.params['N'] * 3.0, places=4)
        # All values should be the same by construction.
        self.assertTrue(all(x - self.res[0, 1] < 1E-8 for x in self.res[:, 1]))


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        TestCylindricalFluxDensityObservable))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
