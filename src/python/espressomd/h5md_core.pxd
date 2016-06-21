from libcpp.string cimport string
cdef extern from "h5md_core.hpp":
    cdef cppclass H5mdCore:
        H5mdCore(const string, const string) except +
        int write_positions()
        int write_velocities()
        int write_forces()

