# Copyright (C) 2019 The ESPResSo project
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

try:
    import MDAnalysis  # pylint: disable=unused-import
except ImportError:
    sample = importlib_wrapper.MagicMock()
    skipIfMissingFeatures = ut.skip(
        "Python module MDAnalysis not available, skipping test!")
else:
    sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
        "@SAMPLES_DIR@/MDAnalysisIntegration.py")


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_universe(self):
        system = self.system
        partcls = system.part.all()
        u = sample.u
        self.assertEqual(len(u.atoms), 10)
        np.testing.assert_equal(u.atoms.ids, np.arange(10) + 1)
        np.testing.assert_equal(u.atoms.types, sorted(5 * ['T0', 'T1']))
        np.testing.assert_almost_equal(
            u.atoms.charges, partcls.q, decimal=6)
        np.testing.assert_almost_equal(
            u.atoms.positions, partcls.pos, decimal=6)
        np.testing.assert_almost_equal(
            u.atoms.velocities, partcls.v, decimal=6)
        np.testing.assert_almost_equal(
            u.atoms.forces, partcls.f, decimal=6)


if __name__ == "__main__":
    ut.main()
