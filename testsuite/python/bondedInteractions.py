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
import espressomd._system as es
from espressomd import code_info
import numpy as np
from espressomd.interactions import FeneBond, HarmonicBond, Dihedral, Angle_Harmonic, Angle_Cosine,\
       Angle_Cossquare, Subt_Lj, Stretching_Force, Area_Force_Local, Bending_Force, Volume_Force,\
       Area_Force_Global, Stretchlin_Force




class ParticleProperties(ut.TestCase):
#  def __init__(self,particleId):
#    self.pid=particleId
  # the system which will be tested
  system = espressomd.System()
  
  # Particle id to work on
  pid=17
  
  # Error tolerance when comparing arrays/tuples...
  tol=1E-9

  def bondsMatch(self,inType,outType,inParams,outParams):
    """Check, if the bond type set and gotten back as well as the bond 
    parameters set and gotten back match. Only check keys present in
    inParams.
    """
    if inType!=outType:
      return False

    for k in inParams.keys():
      if k not in outParams:
        return False
      if outParams[k]!=inParams[k]:
        return False
    
    return True




  def setUp(self):
      self.system.part[self.pid].pos =0,0,0

  def generateTestForBondParams(_bondId,_bondClass,_params):
    """Generates test cases for checking bond parameters set and gotten back
    from Es actually match. Only keys which are present  in _params are checked
    1st arg: Id of the bonded ia in Espresso to test on, i.e., 0,2,1...
    2nd: Class of the bond potential to test, ie.e, FeneBond, HarmonicBond
    3rd: Bond parameters as dictionary, i.e., {"k"=1.,"r_0"=0.
    """
    bondId=_bondId
    bondClass=_bondClass
    params=_params


    def func(self):
        # This code is run at the execution of the generated function.
        # It will use the state of the variables in the outer function, 
        # which was there, when the outer function was called
        self.system.bondedInter[bondId]=bondClass(**params)
        outBond=self.system.bondedInter[bondId]
        tnIn=bondClass(**params).typeNumber()
        tnOut=outBond.typeNumber()
        outParams=outBond.params
        self.assertTrue(self.bondsMatch(tnIn,tnOut,params,outParams), bondClass(**params).typeName()+": value set and value gotten back differ for bond id "+str(bondId)+": "+params.__str__()+" vs. "+outParams.__str__())

    return func
  

  def test_aa_bondedInterSetterGetter(self):
    self.system.bondedInter[0]=HarmonicBond(k=0,r_0=0)
    bond=self.system.bondedInter[0]
    self.assertTrue(isinstance(bond,HarmonicBond),"The bond was created as harmonic bond but the instance gotten back is of different type.")
  
  test_fene = generateTestForBondParams(0,FeneBond,{"r_0":1.1, "k":5.2, "d_r_max":3.})
  test_fene2 = generateTestForBondParams(1,FeneBond,{"r_0":1.1, "k":5.2, "d_r_max":3.})
  test_harmonic = generateTestForBondParams(0,HarmonicBond,{"r_0":1.1, "k":5.2})
  test_harmonic2 = generateTestForBondParams(0,HarmonicBond,{"r_0":1.1, "k":5.2, "r_cut":1.3})
  test_dihedral = generateTestForBondParams(0,Dihedral,{"mult":3.0, "bend":5.2,"phase":3.})
  test_angle_harm = generateTestForBondParams(0,Angle_Harmonic,{"bend":5.2, "phi0":3.2})
  test_angle_cos = generateTestForBondParams(0,Angle_Cosine,{"bend":5.2, "phi0":3.2})
  test_angle_cossquare = generateTestForBondParams(0,Angle_Cossquare,{"bend":5.2, "phi0":0.})
  test_subt_lj = generateTestForBondParams(0,Subt_Lj,{"k":5.2, "r":3.2})
  test_stretching_force = generateTestForBondParams(0,Stretching_Force,{"r0":5.2, "ks":3.2})
  test_area_force_local = generateTestForBondParams(0,Area_Force_Local,{"A0_l":5.2, "ka_l":3.2})
  test_bending_force = generateTestForBondParams(0,Bending_Force,{"phi0":5.2, "kb":3.2})
  test_volume_force = generateTestForBondParams(0,Volume_Force,{"V0":5.2, "kv":3.2})
  test_area_force_global = generateTestForBondParams(0,Area_Force_Global,{"A0_g":5.2, "ka_g":3.2})
  test_stretchlin_force = generateTestForBondParams(0,Stretchlin_Force,{"r0":5.2, "kslin":3.2})


if __name__ == "__main__":
 print("Features: ",code_info.features())
 ut.main()
