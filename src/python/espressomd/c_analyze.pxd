# For C-extern Analysis 

from System cimport *
cimport numpy as np
from utils cimport *

cdef extern from "statistics.hpp":
  ctypedef struct Observable_stat:
    int init_status
    DoubleList data
    int n_coulomb
    int n_dipolar
    int n_non_bonded
    double *bonded
    double *non_bonded
    double *coulomb
    double *dipolar

cdef extern from "statistics.hpp":
  cdef double mindist(IntList *set1, IntList *set2)
  cdef double distto(double pos[3], int pid)
  cdef double *obsstat_bonded(Observable_stat *stat, int j)
  cdef double *obsstat_nonbonded(Observable_stat *stat, int i, int j)

cdef extern from "energy.hpp":
  cdef Observable_stat total_energy

cdef extern from "energy.hpp":
  cdef void master_energy_calc()
  cdef void init_energies(Observable_stat *stat)
