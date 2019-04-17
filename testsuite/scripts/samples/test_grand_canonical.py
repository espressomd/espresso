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
    "@SAMPLES_DIR@/grand_canonical.py", n_int_cycles=51, n_int_steps=5,
    warm_n_times=10, warm_steps=50, cmd_arguments=[0.01, 0.01])


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_deviation_to_target_concentration(self):
        # deviation < 10%
        self.assertLess(abs(sample.deviation), 10.0)


if __name__ == "__main__":
    ut.main()
