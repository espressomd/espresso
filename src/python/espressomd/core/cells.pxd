from libcpp.vector cimport vector
from libcpp.pair cimport pair

cdef extern from "cells.hpp":
    int CELL_STRUCTURE_CURRENT
    int CELL_STRUCTURE_DOMDEC
    int CELL_STRUCTURE_NSQUARE
    int CELL_STRUCTURE_LAYERED

    vector[pair[int, int]] mpi_get_pairs(double distance)
