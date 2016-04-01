#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
import espressomd
import unittest as ut
import numpy as np

if "EWALD_GPU" in espressomd.features():
    from espressomd.electrostatics import EwaldGpu

    class ewald_GPU_test(ut.TestCase):
        es = espressomd.System()
        test_params = {}
        test_params["bjerrum_length"] = 2
        test_params["num_kx"] = 2
        test_params["num_ky"] = 2
        test_params["num_kz"] = 2
        test_params["K_max"] = 10
        test_params["time_calc_steps"] = 100
        test_params["rcut"] = 0.9
        test_params["accuracy"] = 1e-1
        test_params["precision"] = 1e-2
        test_params["alpha"] = 3.5

        def runTest(self):
            ewald = EwaldGpu(bjerrum_length=self.test_params["bjerrum_length"], num_kx=self.test_params["num_kx"], num_ky=self.test_params["num_ky"], num_kz=self.test_params["num_kz"], rcut=self.test_params[
                "rcut"], accuracy=self.test_params["accuracy"], precision=self.test_params["precision"], alpha=self.test_params["alpha"], time_calc_steps=self.test_params["time_calc_steps"], K_max=self.test_params["K_max"])
            self.es.actors.add(ewald)
            set_params = ewald._getParamsFromEsCore()
            SAME = True
            for i in self.test_params.keys():
                if set_params[i] != self.test_params[i]:
                    print "Parameter mismatch: ", i, set_params[i], self.test_params[i]
                    SAME = False
                    break
                    return SAME

    if __name__ == "__main__":
        print("Features: ", espressomd.features())
        ut.main()
