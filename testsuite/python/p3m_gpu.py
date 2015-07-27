#
# Copyright (C) 2013,2014 The ESPResSo project
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
from espressomd.electrostatics import P3M_GPU

class P3M_GPU_test(ut.TestCase):
    es=espressomd.System()
    test_params={}
    test_params["bjerrum_length"]=2
    test_params["cao"]=2
    test_params["inter"]=3
    test_params["r_cut"]=0.9
    test_params["accuracy"]=1e-1
    test_params["mesh"]=[10, 10, 10]
    test_params["epsilon"]=20.0
    test_params["mesh_off"] = [0.8, 0.8, 0.8]
    test_params["alpha"] = 1.1
    test_params["tune"] = False

    def runTest(self):
        p3m=P3M_GPU(bjerrum_length=self.test_params["bjerrum_length"], cao=self.test_params["cao"], inter=self.test_params["inter"], r_cut=self.test_params["r_cut"], accuracy=self.test_params["accuracy"], mesh=self.test_params["mesh"], epsilon=self.test_params["epsilon"], mesh_off=self.test_params["mesh_off"], alpha=self.test_params["alpha"], tune=self.test_params["tune"])
        self.es.Actors.add(p3m)
        set_params=p3m._getParamsFromEsCore()
        SAME=True
        for i in self.test_params.keys():
            if set_params[i] != self.test_params[i]:
                print "Parameter mismatch: ", i
                SAME=False
                break
        return SAME

if __name__ == "__main__":
 print("Features: ",espressomd.features())
 ut.main()
