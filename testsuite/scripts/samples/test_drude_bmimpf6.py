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

gpu = "gpu" in "@TEST_LABELS@".split(";")
if gpu:
    kwargs = {"n_int_steps": 50, "script_suffix": "gpu", "n_ionpairs": 100,
              "cmd_arguments": ["--path", "./bmimpf6_bulk/gpu/", "--gpu"]}
else:
    kwargs = {"n_int_steps": 80, "script_suffix": "cpu", "n_ionpairs": 60,
              "cmd_arguments": ["--path", "./bmimpf6_bulk/cpu/"]}

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/drude_bmimpf6.py", gpu=gpu, n_max_steps=1000,
    n_timing_steps=10, n_int_cycles=5, **kwargs)


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_rdf(self):
        # baseline at large distances should be close to 1.0 +/- 10%
        self.assertLess(abs(np.mean(sample.rdf_00[-20:]) - 1.0), 0.1)
        self.assertLess(abs(np.mean(sample.rdf_01[-20:]) - 1.0), 0.1)
        self.assertLess(abs(np.mean(sample.rdf_11[-20:]) - 1.0), 0.1)


if __name__ == "__main__":
    ut.main()
