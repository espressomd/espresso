# Copyright (C) 2020 The ESPResSo project
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
    "@SAMPLES_DIR@/constant_pH_ensemble_complex_reaction.py", random_seeds=False)


# @skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_concentrations(self):
        type_H = 3
        err_msg = "Concentration of species {} doesn't match analytical result"
        for ptype in sample.types:
            if ptype != type_H:
                self.assertAlmostEqual(sample.concentrations[ptype],
                                       sample.concentrations_numerical[ptype],
                                       delta=1e-3,
                                       msg=err_msg.format(sample.types_name[ptype]))
                self.assertLess(sample.concentrations_95ci[ptype], 1e-3,
                                msg="95% confidence interval too large")
        self.assertAlmostEqual(sample.K_sim, sample.K, delta=1e-3)
        self.assertAlmostEqual(sample.N0_sim, sample.N0, delta=1e-3)


if __name__ == "__main__":
    ut.main()
