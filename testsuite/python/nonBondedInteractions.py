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
from espressomd.interactions import LennardJonesInteraction


class NonBondedInteractionsTests(ut.TestCase):
#  def __init__(self,particleId):
#    self.pid=particleId
  


  # Handle to espresso system
  es=espressomd.System()

  def intersMatch(self,inType,outType,inParams,outParams):
    """Check, if the interaction type set and gotten back as well as the bond 
    parameters set and gotten back match. Only check keys present in
    inParams.
    """
    if inType!=outType:
      print("Type mismatch:",inType,outType)
      return False

    for k in inParams.keys():
      if k not in outParams:
        print(k,"missing from returned parameters")
        return False
      if outParams[k]!=inParams[k]:
        print("Mismatch in parameter ",k,inParams[k],outParams[k])
	return False
    
    return True





  def generateTestForNonBondedInteraction(_partType1,_partType2,_interClass,_params,_interName):
    """Generates test cases for checking interaction parameters set and gotten back
    from Es actually match. Only keys which are present  in _params are checked
    1st and 2nd arg: Particle type ids to check on
    3rd: Class of the interaction to test, ie.e, FeneBond, HarmonicBond
    4th: Interaction parameters as dictionary, i.e., {"k"=1.,"r_0"=0.
    5th: Name of the interaction property to set (i.e. "lennardJones")
    """
    partType1=_partType1
    partType2=_partType2
    interClass=_interClass
    params=_params
    interName=_interName


    def func(self):
      # This code is run at the execution of the generated function.
      # It will use the state of the variables in the outer function, 
      # which was there, when the outer function was called
      
      # Set parameters
      getattr(self.es.nonBondedInter[partType1,partType2],interName).setParams(**params)
      
      # Read them out again
      outInter=getattr(self.es.nonBondedInter[partType1,partType2],interName)
      outParams=outInter.getParams()
      
      
      self.assertTrue(self.intersMatch(interClass,type(outInter),params,outParams), interClass(**params).typeName()+": value set and value gotten back differ for particle types "+str(partType1)+" and "+str(partType2)+": "+params.__str__()+" vs. "+outParams.__str__())

    return func

  test_lj1=generateTestForNonBondedInteraction(\
    0,0,LennardJonesInteraction,\
    {"epsilon":1.,"sigma":2.,"cutoff":3.,"shift":4.,"offset":5.,"min":7.},\
    "lennardJones")
  test_lj2=generateTestForNonBondedInteraction(\
    0,0,LennardJonesInteraction,\
    {"epsilon":1.3,"sigma":2.2,"cutoff":3.4,"shift":4.1,"offset":5.1,"min":7.1},\
    "lennardJones")
  test_lj3=generateTestForNonBondedInteraction(\
    0,0,LennardJonesInteraction,\
    {"epsilon":1.3,"sigma":2.2,"cutoff":3.4,"shift":4.1,"offset":5.1,"min":7.1},\
    "lennardJones")
    
  test_ljgen1=generateTestForNonBondedInteraction(\
    0,0,GenericLennardJonesInteraction,\
    {"epsilon":1.,"sigma":2.,"cutoff":3.,"shift":4.,"offset":5.,"e1":7,"e2":8,"b1":9.,"b2":10.,"lambda":11.,"delta":12.},\
    "genericLennardJones")
  test_ljgen2=generateTestForNonBondedInteraction(\
    0,0,GenericLennardJonesInteraction,\
    {"epsilon":1.1,"sigma":2.1,"cutoff":3.1,"shift":4.1,"offset":5.1,"e1":71,"e2":81,"b1":9.1,"b2":10.1,"lambda":11.1,"delta":12.1},\
    "genericLennardJones")
  test_ljgen3=generateTestForNonBondedInteraction(\
    0,0,GenericLennardJonesInteraction,\
    {"epsilon":1.2,"sigma":2.2,"cutoff":3.2,"shift":4.2,"offset":5.2,"e1":72,"e2":82,"b1":9.2,"b2":10.2,"lambda":11.2,"delta":12.2},\
    "genericLennardJones")

  def test_forcecap(self):
    self.es.nonBondedInter.setForceCap(17.5)
    self.assertEqual(self.es.nonBondedInter.getForceCap(),17.5)


if __name__ == "__main__":
 print("Features: ",espressomd.features())
 ut.main()

