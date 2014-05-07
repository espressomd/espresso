include "myconfig.pxi"
cimport cuda_init
    
cdef class CudaInitHandle:
  def __init__(self):
    IF CUDA != 1:
      raise Exception("Cuda is not compiled in")
  
  property device:
    IF CUDA == 1:
      def __set__(self, int _dev):
        if cuda_set_device(_dev):
          raise Exception("cuda device set error")
      def __get__(self): 
        dev = cuda_get_device()
        if dev == -1:
          raise Exception("cuda device get error")
        return dev

  # property device_list:
  #   IF CUDA == 1:
  #     def __set__(self, int _dev):
  #       raise Exception("cuda device list is read only")
  #     def __get__(self): 
  #       cdef int _p_devl
  #       cdef char _devname[4+64]
  #       if getdevicelist(&_p_devl, _devname):
  #         raise Exception("cuda devicelist error")
  #       return _devname


