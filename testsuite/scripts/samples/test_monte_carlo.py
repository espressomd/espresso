#
# Copyright (C) 2023 The ESPResSo project
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
import importlib_wrapper

mode, method = "@TEST_SUFFIX@".split("_")
assert method in ("cph", "re")

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/monte_carlo.py", script_suffix="@TEST_SUFFIX@",
    cmd_arguments=["--mode", mode, "--method", method], sample_size=150)


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test(self):
        if method == "cph":
            self.assertAlmostEqual(sample.alpha_avg, 0.29, delta=0.05)
            self.assertAlmostEqual(sample.alpha_err, 0.01, delta=0.02)
            self.assertAlmostEqual(sample.acceptance_rate, 0.56, delta=0.10)
        else:
            self.assertAlmostEqual(sample.alpha_avg, 0.33, delta=0.05)
            self.assertAlmostEqual(sample.alpha_err, 0.01, delta=0.02)
            self.assertAlmostEqual(sample.acceptance_rate, 0.55, delta=0.10)


if __name__ == "__main__":
    ut.main()
