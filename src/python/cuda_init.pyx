include "myconfig.pxi"
cimport cuda_init

    
cdef class CudaInitHandle:
  def __init__(self):
    IF LB_GPU != 1:
      raise Exception("Cuda is not compiled in")
  
  property device:
    IF LB_GPU == 1:
      def __set__(self, int _dev):
        if setdevice(_dev):
          raise Exception("cuda device set error")
      def __get__(self): 
        cdef int _p_dev
        if getdevice(&_p_dev):
          raise Exception("cuda device get error")
        return _p_dev

  property device_list:
    IF LB_GPU == 1:
      def __set__(self, int _dev):
        raise Exception("cuda device list is read only")
      def __get__(self): 
        cdef int _p_devl
        cdef char _devname[4+64]
        if getdevicelist(&_p_devl, _devname):
          raise Exception("cuda devicelist error")
        return _devname


