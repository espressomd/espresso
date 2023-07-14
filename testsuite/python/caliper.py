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
import unittest_decorators as utx
import espressomd.profiler
import subprocess
import pathlib
import sys
import os

EXPECTED_LABELS = """
calc_energies
  short_range_loop
integrate
  Integration loop
    force_calc
      copy_forces_from_GPU
      short_range_loop
      calc_long_range_forces
      init_forces
      copy_particles_to_GPU
"""


@utx.skipIfMissingFeatures(["CALIPER"])
class Test(ut.TestCase):

    @utx.skipIfMissingFeatures(["P3M", "WCA"])
    def test_runtime_report(self):
        has_cuda = espressomd.has_features(["CUDA"])
        script = str(pathlib.Path(__file__).parent / "caliper_child.py")
        my_env = os.environ.copy()
        my_env["CALI_CONFIG_PROFILE"] = "runtime-report"
        process = subprocess.Popen([sys.executable, script],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   env=my_env)
        stderr = process.communicate(timeout=60)[1].decode("utf-8")
        lines = stderr.strip().split("\n")
        for i, line in enumerate(lines):
            if not line.startswith("WARNING:"):
                lines = lines[i:]
                break
        header = "Path\tInclusive time\tExclusive\ttime\tTime %"
        self.assertEqual(lines[0].split(), header.split(),
                         msg=f"Caliper summary should start with '{header}'")
        labels = [line[:30].strip() for line in lines[1:]]
        labels_ref = [x.strip() for x in EXPECTED_LABELS.strip().split("\n")
                      if "GPU" not in x.upper() or has_cuda]
        self.assertEqual(labels[:len(labels_ref)], labels_ref,
                         msg=f"Caliper returned this summary:\n{stderr}")

    def test_interface(self):
        cali = espressomd.profiler.Caliper()
        self.assertIsNone(cali.call_method("unknown"))


if __name__ == "__main__":
    ut.main()
