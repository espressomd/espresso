from __future__ import print_function, absolute_import

cdef extern from "communication.hpp":

    void mpi_lees_edwards_image_reset()
