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
import unittest as ut
import espressomd
import numpy as np
from espressomd.magnetostatics import *


class MagnetostaticsInteractionsTests(ut.TestCase):
    # Handle to espresso system
    system=espressomd.System()

    def paramsMatch(self,inParams,outParams):
        """Check, if the parameters set and gotten back match.
        Only check keys present in inParams.
        """

        # Delete the bjerrum_length from outparams, because
	# it is calculated from the prefactor in some cases
	# and therefore not necessarily in inparams.
	if "bjerrum_length" in outParams:
	    del outParams["bjerrum"]

	for k in inParams.keys():
            if k not in outParams:
                print(k,"missing from returned parameters")
                return False
            if outParams[k]!=inParams[k]:
                print("Mismatch in parameter ",k,inParams[k],outParams[k])
                return False
    
        return True


    def setUp(self):
        self.system.box_l = 10,10,10
        self.system.part[0].pos =0,0,0
        self.system.part[1].pos =0.1,0.1,0.1
        self.system.part[0].dip =1.3,2.1,-6
        self.system.part[1].dip =4.3,4.1,7.5
      

    def generateTestForMagnetostaticInteraction(_interClass,_params):
        """Generates test cases for checking interaction parameters set and gotten back
        from Es actually match. Only keys which are present  in _params are checked
        1st: Interaction parameters as dictionary, i.e., {"k"=1.,"r_0"=0.
        2nd: Name of the interaction property to set (i.e. "P3M")
        """
        params=_params
        interClass=_interClass


        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function, 
            # which was there, when the outer function was called
            
            # set Parameter
            Inter=interClass(**params)
            Inter.validateParams()
            Inter._setParamsInEsCore()
            
            # Read them out again
            outParams=Inter.getParams()
      
      
            self.assertTrue(self.paramsMatch(params,outParams), "Missmatch of parameters.\nParameters set "+params.__str__()+" vs. output parameters "+outParams.__str__())

        return func

    test_DP3M=generateTestForMagnetostaticInteraction(DipolarP3M,dict(prefactor = 1.0,\
                                                              epsilon = 0.0,\
                                                              inter = 1000,\
                                                              mesh_off = [0.5,0.5,0.5],\
                                                              r_cut = 2.4,\
                                                              mesh = [8,8,8],\
                                                              cao = 1,\
                                                              alpha = 12,\
                                                              accuracy = 0.01))


    test_DdsCpu=generateTestForMagnetostaticInteraction(DipolarDirectSumCpu,dict(prefactor = 3.4))
    test_DdsRCpu=generateTestForMagnetostaticInteraction(DipolarDirectSumWithReplicaCpu,dict(prefactor = 3.4,n_replica=2))


if __name__ == "__main__":
    print("Features: ",espressomd.features())
    ut.main()

