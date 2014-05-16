# Tests particle property setters/getters
import unittest as ut
import espresso.System as es
import numpy as np
from espresso.interactions import FeneBond



class ParticleProperties(ut.TestCase):
#  def __init__(self,particleId):
#    self.pid=particleId
  
  # Particle id to work on
  pid=17
  
  # Error tolerance when comparing arrays/tuples...
  tol=1E-9

  def arraysNearlyEqual(self,a,b):
    """Test, if the magnitude of the difference between two arrays is smaller than the tolerance"""

    # Check length
    if len(a) != len(b):
      return False


    # We have to use a loop, since we can't be sure, we're getting numpy arrays
    sum=0.
    for i in range(len(a)):
      sum+= abs(a[i]-b[i])

    if sum >self.tol:
      return False

    return True
      



  def setUp(self):
    es.part[self.pid].pos =0,0,0
    es.bondedInter[0]=FeneBond(k=1,d_r_max=5)
    es.bondedInter[1]=FeneBond(k=1,d_r_max=5)

  def generateTestForVectorProperty(_propName,_value):
    """Generates test cases for vectorial particle properties such as
    position, velocity...
    1st arg: name of the property (e.g., "pos"), 
    2nd array: value to be used for testing. Has to be numpy.array of floats
    """
    # This is executed, when generateTestForVectorProperty() is called
    propName=_propName
    value=_value

    def func(self):
      # This code is run at the execution of the generated function.
      # It will use the state of the variables in the outer function, 
      # which was there, when the outer function was called
      setattr(es.part[self.pid],propName,value)
      print(propName,value,getattr(es.part[self.pid],propName))
      self.assertTrue(self.arraysNearlyEqual(getattr(es.part[self.pid],propName), value),propName+": value set and value gotten back differ.")

    return func
  
  def generateTestForScalarProperty(_propName,_value):
    """Generates test cases for scalar particle properties such as
    type, mass, charge...
    1st arg: name of the property (e.g., "type"), 
    2nd array: value to be used for testing. int or float
    """
    # This is executed, when generateTestForVectorProperty() is called
    propName=_propName
    value=_value

    def func(self):
      # This code is run at the execution of the generated function.
      # It will use the state of the variables in the outer function, 
      # which was there, when the outer function was called
      setattr(es.part[self.pid],propName,value)
      print(propName,value,getattr(es.part[self.pid],propName))
      self.assertTrue(getattr(es.part[self.pid],propName)==value,propName+": value set and value gotten back differ.")

    return func


  test_pos=generateTestForVectorProperty("pos",np.array([0.1,0.2,0.3]))
  test_v=generateTestForVectorProperty("v",np.array([0.2,0.3,0.4]))
  test_f=generateTestForVectorProperty("f",np.array([0.2,0.3,0.7]))
  test_type=generateTestForScalarProperty("type",int(3))
  test_bonds_property=generateTestForScalarProperty("bonds", ((0,1),(1,2)))


  if "MASS" in es.code_info.features(): 
    test_mass=generateTestForScalarProperty("mass",1.3)

  if "ROTATION" in es.code_info.features(): 
    test_omega_lab=generateTestForVectorProperty("omega_lab",np.array([4.,2.,1.]))
    test_omega_body=generateTestForVectorProperty("omega_body",np.array([4.,72.,1.]))
    test_torque_lab=generateTestForVectorProperty("torque_lab",np.array([4.,72.,3.7]))
    # The tested value has to be nromalized!
    test_quat=generateTestForVectorProperty("quat",np.array([0.5,0.5,0.5,0.5]))
#    test_director=generateTestForVectorProperty("director",np.array([0.5,0.4,0.3]))


  if "ELECTROSTATICS" in es.code_info.features():
    test_charge=generateTestForScalarProperty("q",-19.7)

  if "DIPOLES" in es.code_info.features():
    test_dip=generateTestForVectorProperty("dip",np.array([0.5,-0.5,3]))
    test_dipm=generateTestForScalarProperty("dipm",-9.7)

  if "VIRTUAL_SITES" in es.code_info.features():
    test_virtual=generateTestForScalarProperty("virtual",1)
  if "VIRTUAL_SITES_RELATIVE" in es.code_info.features():
    test_zz_vs_relative=generateTestForScalarProperty("vs_relative",((0,5.0)))
    

if __name__ == "__main__":
 print("Features: ",es.code_info.features())
 ut.main()
