#
# Copyright (C) 2021-2022 The ESPResSo project
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
#
import unittest as ut
import espressomd
import numpy as np


class BoxGeometry(ut.TestCase):
    box_l = [5.0, 5.0, 5.0]
    system = espressomd.System(box_l=box_l)

    def test_box_length_interface(self):
        for i in range(3):
            local_box_l = self.box_l.copy()
            local_box_l[i] = -1.

            with self.assertRaisesRegex(ValueError, "Attribute 'box_l' must be > 0"):
                self.system.box_l = local_box_l

            # the box length should not be updated
            np.testing.assert_equal(self.box_l, np.copy(self.system.box_l))

        with self.assertRaisesRegex(RuntimeError, "Provided argument of type 'std::vector<double>{.size=2}' is not convertible to 'Utils::Vector<double, 3>'"):
            self.system.box_l = self.box_l[:2]

    def test_periodicity(self):
        import itertools
        for periodicity in itertools.product((True, False), repeat=3):
            self.system.periodicity = periodicity

            np.testing.assert_equal(np.copy(self.system.periodicity),
                                    periodicity)

        default_periodicity = (True, True, True)
        self.system.periodicity = default_periodicity

        with self.assertRaisesRegex(RuntimeError, "Provided argument of type 'std::vector<ScriptInterface::Variant{bool}>{.size=2}' is not convertible to 'Utils::Vector<bool, 3>'"):
            self.system.periodicity = (True, True)

        # the periodicity should not be updated
        np.testing.assert_equal(np.copy(self.system.periodicity),
                                default_periodicity)

    def test_min_global_cut(self):
        self.system.min_global_cut = 0.01
        with self.assertRaisesRegex(ValueError, "Attribute 'min_global_cut' must be >= 0"):
            self.system.min_global_cut = -2.

        # the min global cutoff should not be updated
        np.testing.assert_equal(self.system.min_global_cut, 0.01)


if __name__ == "__main__":
    ut.main()
