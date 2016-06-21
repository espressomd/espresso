cdef class PyH5mdCore:
    cdef H5mdCore* c_h5md
    def __cinit__(self, filename, path_to_pythonscript):
        self.c_h5md = new H5mdCore(filename, path_to_pythonscript)
    def __dealloc__(self):
        del self.c_h5md
