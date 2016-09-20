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
from __future__ import print_function
import espressomd
import unittest as ut
import numpy as np
from tests_common import *

if "ELECTROSTATICS" in espressomd.features() and "CUDA" in espressomd.features():
    from espressomd.electrostatics import P3M_GPU

    class P3M_GPU_test(ut.TestCase):
        def runTest(self):
            es = espressomd.System()
            test_params = {}
            test_params["bjerrum_length"] = 2
            test_params["cao"] = 2
            test_params["inter"] = 3
            test_params["r_cut"] = 0.9
            test_params["accuracy"] = 1e-1
            test_params["mesh"] = [10, 10, 10]
            test_params["epsilon"] = 20.0
            test_params["mesh_off"] = [0.8, 0.8, 0.8]
            test_params["alpha"] = 1.1
            test_params["tune"] = False
    
            p3m = P3M_GPU(**test_params)
            es.actors.add(p3m)
            self.assertTrue(params_match(test_params,p3m._get_params_from_es_core()))

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
