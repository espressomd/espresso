# Copyright (C) 2019-2022 The ESPResSo project
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

import unittest as ut
import importlib_wrapper
import numpy as np

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/lennard_jones/lennard_jones.py", N_SAMPLES=300)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):

    def test_rdf(self):
        np.testing.assert_allclose(tutorial.rdf, tutorial.theo_rdf,
                                   rtol=0.0, atol=0.1)

    def test_potential_energy(self):
        # Test that the potential energy/particle agrees with
        # the value from Verlet, Phys. Rev. 1967
        ref_energy = -5.38
        self.assertAlmostEqual(
            tutorial.mean_pot_energy_corrected, ref_energy, 1)


if __name__ == "__main__":
    ut.main()
