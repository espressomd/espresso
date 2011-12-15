
cimport global_variables
import numpy as np

cpdef int marcello_test(int x):
  return x**2

cdef class GlobalsHandle:
  def __init__(self):
    pass

  property time_step:
    def __set__(self, double _time_step):
      mpi_set_time_step(_time_step)
    def __get__(self):
      global time_step
      return time_step

  property skin:
    def __set__(self, double _skin):
      global skin
      skin=_skin
      mpi_bcast_parameter(28)
    def __get__(self): 
      global skin
      return skin

  property box_l:
    def __set__(self, _box_l):
      global box_l
      for i in range(3):
        box_l[i]=_box_l[i]
      mpi_bcast_parameter(0)
    def __get__(self):
      return np.array([box_l[0], box_l[1], box_l[2]])




