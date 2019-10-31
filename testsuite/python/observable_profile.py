# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import sys
import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd
import espressomd.observables


class ProfileObservablesTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.box_l = [10.0, 10.0, 10.0]
    system.cell_system.skin = 0.1
    system.time_step = 0.01
    system.part.add(id=0, pos=[4.0, 4.0, 6.0], v=[0.0, 0.0, 1.0])
    system.part.add(id=1, pos=[4.0, 4.0, 6.0], v=[0.0, 0.0, 1.0])
    flat_index = np.ravel_multi_index((0, 0, 1), (2, 2, 2))
    flat_index_3d = np.ravel_multi_index((0, 0, 1, 2), (2, 2, 2, 3))
    bin_volume = 5.0**3
    kwargs = {'ids': [0, 1],
              'n_x_bins': 2,
              'n_y_bins': 2,
              'n_z_bins': 2,
              'min_x': 0.0,
              'max_x': 10.0,
              'min_y': 0.0,
              'max_y': 10.0,
              'min_z': 0.0,
              'max_z': 10.0}

    def test_density_profile(self):
        density_profile = espressomd.observables.DensityProfile(**self.kwargs)
        self.assertEqual(density_profile.calculate()[
                         self.flat_index], 2.0 / self.bin_volume)
        self.assertEqual(density_profile.n_values(),
                         self.kwargs['n_x_bins'] * self.kwargs['n_y_bins'] * self.kwargs['n_z_bins'])

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_force_density_profile(self):
        density_profile = espressomd.observables.ForceDensityProfile(
            **self.kwargs)
        self.system.part[0].ext_force = [0.0, 0.0, 1.0]
        self.system.part[1].ext_force = [0.0, 0.0, 1.0]
        self.system.integrator.run(0)
        self.assertEqual(density_profile.calculate()[
                         self.flat_index_3d], 2.0 / self.bin_volume)
        self.assertEqual(density_profile.n_values(),
                         3 * self.kwargs['n_x_bins'] * self.kwargs['n_y_bins'] * self.kwargs['n_z_bins'])

    def test_flux_density_profile(self):
        density_profile = espressomd.observables.FluxDensityProfile(
            **self.kwargs)
        self.system.integrator.run(0)
        self.assertEqual(density_profile.calculate()[
                         self.flat_index_3d], 2.0 / self.bin_volume)
        self.assertEqual(density_profile.n_values(),
                         3 * self.kwargs['n_x_bins'] * self.kwargs['n_y_bins'] * self.kwargs['n_z_bins'])


if __name__ == '__main__':
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        ProfileObservablesTest))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
