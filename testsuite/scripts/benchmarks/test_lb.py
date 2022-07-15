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

# make simulation deterministic
np.random.seed(42)

benchmark, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@BENCHMARKS_DIR@/lb.py", cmd_arguments=["--particles_per_core", "80"],
    measurement_steps=200, n_iterations=2)


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = benchmark.system


if __name__ == "__main__":
    ut.main()
