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
    "@SAMPLES_DIR@/gibbs_ensemble/gibbs_ensemble_socket.py",
    cmd_arguments=["-S", "2000", "-C", "@SAMPLES_DIR@/gibbs_ensemble/lj_gibbs_client.py",
                   "-E", "@CMAKE_BINARY_DIR@/pypresso", "-T", "0.1", "-s", "1"])


@skipIfMissingFeatures
class Sample(ut.TestCase):
    sample

    def test(self):
        self.assertEqual(sample.kT, 0.1) # specified in cmake list
        
        # Check densities. The system will not have converged after 1k steps,
        # but the two densities should be subtantially different
        last_densities = [sample.densities[0][-1], sample.densities[1][-1]]
        self.assertLess(min(last_densities), 0.25)
        self.assertGreater(max(last_densities), 0.4)



if __name__ == "__main__":
    ut.main()
