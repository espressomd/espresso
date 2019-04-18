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

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/diffusion_coefficient.py")


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_Green_Kubo_accuracy(self):
        # test simulated value to +/- 15% accuracy
        self.assertGreater(sample.ratio, 0.85)
        self.assertLess(sample.ratio, 1.15)


if __name__ == "__main__":
    ut.main()
