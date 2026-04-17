#
# Copyright (C) 2023-2026 The ESPResSo project
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
integrate
  update_cabana_state
  Initial Force Calculation
    calculate_forces
      copy_particles_to_GPU
      update_cabana_state
      init_forces_and_thermostat
      calc_long_range_forces
      cabana_short_range
      copy_forces_from_GPU
  Integration loop
    integrator_step_1
    resort_particles_if_needed
    calculate_forces
      copy_particles_to_GPU
      update_cabana_state
      init_forces_and_thermostat
      calc_long_range_forces
      cabana_short_range
      copy_forces_from_GPU
    integrator_step_2
calc_energies
  update_cabana_state
  cabana_short_range
"""


@utx.skipIfMissingFeatures(["CALIPER"])
class Test(ut.TestCase):

    @utx.skipIfMissingFeatures(["P3M", "WCA"])
    def test_runtime_report(self):
        has_cuda = espressomd.has_features(["CUDA"])
        script = str(pathlib.Path(__file__).parent / "caliper_child.py")
        my_env = os.environ.copy()
        my_env["CALI_CONFIG"] = "runtime-report"
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
        header = "Path\tMin time/rank\tMax time/rank\tAvg time/rank\tTime %"
        self.assertEqual(lines[0].split(), header.split(),
                         msg=f"Caliper summary should start with '{header}'")
        labels = [line[:36].rstrip() for line in lines[1:]]
        labels_ref = [x.rstrip() for x in EXPECTED_LABELS.strip().split("\n")
                      if x.rstrip() and ("GPU" not in x.upper() or has_cuda)]
        for l in labels_ref: assert l in labels, f"Label {l} not in profile"

    def test_interface(self):
        cali = espressomd.profiler.Caliper()
        self.assertIsNone(cali.call_method("unknown"))


if __name__ == "__main__":
    ut.main()
