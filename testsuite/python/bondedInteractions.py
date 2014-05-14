# Tests particle property setters/getters
import unittest as ut
import espresso.System as es
import numpy as np
from espresso.interactions import FeneBond, HarmonicBond



class ParticleProperties(ut.TestCase):
#  def __init__(self,particleId):
#    self.pid=particleId
  
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
    es.part[self.pid].pos =0,0,0

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
      es.bondedInter[bondId]=bondClass(**params)
      outBond=es.bondedInter[bondId]
      tnIn=bondClass(**params).typeNumber()
      tnOut=outBond.typeNumber()
      outParams=outBond.params
      self.assertTrue(self.bondsMatch(tnIn,tnOut,params,outParams), bondClass(**params).typeName()+": value set and value gotten back differ for bond id "+str(bondId)+": "+params.__str__()+" vs. "+outParams.__str__())

    return func
  

  def test_aa_bondedInterSetterGetter(self):
    es.bondedInter[0]=HarmonicBond(k=0,r_0=0)
    bond=es.bondedInter[0]
    self.assertTrue(isinstance(bond,HarmonicBond),"The bond was created as harmonic bond but the instance gotten back is of different type.")
  
  test_harmonic=generateTestForBondParams(0,HarmonicBond,{"r_0":1.1, "k":5.2})
  test_harmonic2=generateTestForBondParams(0,HarmonicBond,{"r_0":1.1, "k":5.2,"r_cut":1.3})
  test_fene=generateTestForBondParams(0,FeneBond,{"r_0":1.1, "k":5.2,"d_r_max":3.})
  test_fene2=generateTestForBondParams(1,FeneBond,{"r_0":1.1, "k":5.2,"d_r_max":3.})


if __name__ == "__main__":
 print("Features: ",es.code_info.features())
 ut.main()
