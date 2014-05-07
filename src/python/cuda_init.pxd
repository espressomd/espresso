cdef extern from "cuda_init.hpp":
    int cuda_set_device(int dev)
    int cuda_get_device()
#    int getdevicelist(int* devl, char* devname)
