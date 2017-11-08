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
        'center': [5.0, 5.0, 0.0],
        'n_r_bins': 100,
        'n_phi_bins': 1,
        'n_z_bins': 1,
        'min_r': 0.0,
        'min_phi': -np.pi,
        'min_z': 0.0,
        'max_r': 5.0,
        'max_phi': np.pi,
        'max_z': 10.0,
        'N': 100
    }

    def get_ind(self, ind):
        return utils.get_unravelled_index(
            [self.params['n_r_bins'], self.params['n_phi_bins'], self.params['n_z_bins'], 3], 4, ind)

    def test_normalization(self):
        self.system.part.clear()
        for i in range(self.params['N']):
            self.system.part.add(id=i, pos=np.array(
                self.params['center']) + np.array([i * 5.0 / self.params['N'], 0, 0]), v=[0.0, 0.0, 3.0])
        # Set up the Observable.
        p = espressomd.observables.CylindricalFluxDensityProfile(
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

        r_bin_size = (
            self.params['max_r'] - self.params['min_r']) / self.params['n_r_bins']
        hist = np.array(p.calculate())
        res = []
        for ind, v_dens in enumerate(hist):
            # Get the unflattened index.
            _ind = self.get_ind(ind)
            # If we have the 'z' component add the velocity density to the
            # result.
            if _ind[3] == 2:
                # Append the mid of the bin position and the z-component of the
                # velocity density.
                res.append([(_ind[0] + 0.5) * r_bin_size, v_dens])

        res = np.array(res)
        # The Riemann sum should yield the sum of all z-velocities
        riemann_sum = np.sum(res[:, 0] * res[:, 1]) * 2.0 * np.pi * \
            (self.params['max_z'] - self.params['min_z']) * r_bin_size
        self.assertAlmostEqual(riemann_sum, self.params['N'] * 3.0, places=4)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        TestCylindricalFluxDensityObservable))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
