import sys
import unittest as ut
import numpy as np
import espressomd
import espressomd.observables


class ProfileObservablesTest(ut.TestCase):
    system = espressomd.System()
    system.box_l = [10.0, 10.0, 10.0]
    system.cell_system.skin = 0.1
    system.time_step = 0.01
    system.part.add(
        id=0, pos=[
            4.0, 4.0, 6.0], v=[
            0.0, 0.0, 1.0], ext_force=[
                0.0, 0.0, 1.0])
    system.part.add(
        id=1, pos=[
            4.0, 4.0, 6.0], v=[
            0.0, 0.0, 1.0], ext_force=[
                0.0, 0.0, 1.0])
    flat_index = np.ravel_multi_index((0, 0, 1), (2, 2, 2))
    flat_index_3d = np.ravel_multi_index((0, 0, 1, 2), (2, 2, 2, 3))
    bin_volume = 5.0**3
    kwargs = {'ids': [0, 1],
              'xbins': 2,
              'ybins': 2,
              'zbins': 2,
              'minx': 0.0,
              'maxx': 10.0,
              'miny': 0.0,
              'maxy': 10.0,
              'minz': 0.0,
              'maxz': 10.0}

    def test_density_profile(self):
        density_profile = espressomd.observables.DensityProfile(**self.kwargs)
        self.assertEqual(density_profile.calculate()[
                         self.flat_index], 2.0 / (2.0 * self.bin_volume))

    def test_force_density_profile(self):
        density_profile = espressomd.observables.ForceDensityProfile(
            **self.kwargs)
        self.system.integrator.run(0)
        self.assertEqual(density_profile.calculate()[
                         self.flat_index_3d], 2.0 / (2.0 * self.bin_volume))

    def test_flux_density_profile(self):
        density_profile = espressomd.observables.FluxDensityProfile(
            **self.kwargs)
        self.system.integrator.run(0)
        self.assertEqual(density_profile.calculate()[
                         self.flat_index_3d], 2.0 / (2.0 * self.bin_volume))

if __name__ == '__main__':
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        ProfileObservablesTest))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
    ut.main()
