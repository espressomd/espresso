from __future__ import print_function
import ctypes
import sys

sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

# Tests particle property setters/getters
import unittest as ut
import espresso as es
import numpy as np



class ParticleProperties(ut.TestCase):
#  def __init__(self,particleId):
#    self.pid=particleId
  pid=17

#  def setUp(self):
#    es.part[self.pid].pos =0,0,0

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
      self.assertTrue(np.all(getattr(es.part[self.pid],propName)==value),propName+": value set and value gotten back differ.")

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
      self.assertTrue(getattr(es.part[self.pid],propName)==value,propName+": value set and value gotten back differ.")

    return func


     
  
  
  
  test_pos=generateTestForVectorProperty("pos",np.array([0.1,0.2,0.3]))
  test_v=generateTestForVectorProperty("v",np.array([0.1,0.2,0.3]))
  test_type=generateTestForScalarProperty("type",1)


if __name__ == "__main__":
 ut.main()

