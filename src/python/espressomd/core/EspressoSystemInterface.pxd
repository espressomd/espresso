from libcpp cimport bool
from core.SystemInterface cimport SystemInterface

cdef extern from "EspressoSystemInterface.hpp":
    cdef cppclass EspressoSystemInterface(SystemInterface):
        @staticmethod
        EspressoSystemInterface * _Instance()
        bool requestRGpu()
        void update()
