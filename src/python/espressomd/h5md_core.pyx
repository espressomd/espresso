cdef class PyH5mdCore:
    cdef File* c_h5md
    def __cinit__(self, filename, path_to_pythonscript):
        self.c_h5md = new File(filename, path_to_pythonscript)
    def __dealloc__(self):
        del self.c_h5md
