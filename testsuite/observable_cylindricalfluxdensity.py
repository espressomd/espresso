import sys
import numpy as np
import unittest as ut
import espressomd
import espressomd.observables
from espressomd import utils


class TestCylindricalFluxDensityObservable(ut.TestCase):
    system = espressomd.System()
    system.time_step = 0.01
    system.box_l = [10.0, 10.0, 10.0]
    system.cell_system.skin = 0.4

    params = {
        'center': [5.0, 5.0, 0.0], # center of the histogram
        'n_r_bins': 10, # number of bins in r
        'n_phi_bins': 1, # -*- in phi
        'n_z_bins': 1, # -*- in z
        'min_r': 0.0,
        'min_phi': -np.pi,
        'min_z': 0.0,
        'max_r': 5.0,
        'max_phi': np.pi,
        'max_z': 10.0,
        'N': 100 # number of particles
    }

    def get_ind(self, ind):
        return utils.get_unravelled_index(
            [self.params['n_r_bins'], self.params['n_phi_bins'], self.params['n_z_bins'], 3], 4, ind)

    def setUp(self):
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

    def test_normalization(self):
        self.res = []
        for ind, v_dens in enumerate(self.hist):
            # Get the unflattened index.
            _ind = self.get_ind(ind)
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
        nphist, _ = np.histogram(self.system.part[:].pos[:,0] - self.params['center'][0], weights=self.system.part[:].v[:,2], bins=self.params['n_r_bins'], range=(self.params['min_r'], self.params['max_r']))
        self.assertAlmostEqual(riemann_sum, np.sum(nphist), places=4)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        TestCylindricalFluxDensityObservable))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
