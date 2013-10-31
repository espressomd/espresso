# Handling of interactions

from espresso cimport *
cimport numpy as np
from utils cimport *

cdef extern from "interaction_data.hpp":
  ctypedef struct IA_parameters:
    double LJ_eps
    double LJ_sig
    double LJ_cut
    double LJ_shift
    double LJ_offset
    double LJ_capradius
    double LJ_min

  cdef IA_parameters *get_ia_param(int i, int j)

cdef extern from "lj.hpp":
  cdef int lennard_jones_set_params(int part_type_a, int part_type_b,
                                        double eps, double sig, double cut,
                                        double shift, double offset,
                                        double cap_radius, double min)
  cdef int ljforcecap_set_params(double ljforcecap)
