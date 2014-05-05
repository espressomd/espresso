cdef extern from "config.hpp":
    pass

cdef extern from "cuda_init.h":
    int setdevice(int dev)
    int getdevice(int* dev)
    int getdevicelist(int* devl, char* devname)
