cimport cuda_init

    
cdef class CudaInitHandle:
  def __init__(self):
    pass
  
  property device:
    def __set__(self, int _dev):
      if setdevice(_dev):
        raise Exception("cuda device set error")
    def __get__(self): 
      cdef int _p_dev
      if getdevice(&_p_dev):
        raise Exception("cuda device get error")
      return _p_dev
