cimport global_variables
import numpy as np

cpdef int marcello_test(int x):
  return x**2

cdef class GlobalsHandle:
  def __init__(self):
    pass

  property time_step:
    def __set__(self, double _time_step):
      if _time_step <= 0:
        raise ValueError("Time Step must be positive")
      mpi_set_time_step(_time_step)
    def __get__(self):
      global time_step
      return time_step

  property skin:
    def __set__(self, double _skin):
      if _skin <= 0:
        raise ValueError("Skin must be >= 0")
      global skin
      skin=_skin
      mpi_bcast_parameter(28)
    def __get__(self): 
      global skin
      return skin

  property box_l:
    def __set__(self, _box_l):
      global box_l
      if len(_box_l) != 3:
        raise ValueError("Box length must be of length 3")
      for i in range(3):
        if _box_l[i] <= 0:
          raise ValueError("Box length must be > 0 in all directions")
        box_l[i]=_box_l[i]

      mpi_bcast_parameter(0)

    def __get__(self):
      return np.array([box_l[0], box_l[1], box_l[2]])




