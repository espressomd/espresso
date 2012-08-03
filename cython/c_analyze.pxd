# For C-extern Analysis 

from espresso cimport *
cimport numpy as np
from utils cimport *

cdef extern from "../src/statistics.h":
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

cdef extern from "../src/statistics.h":
	cdef double mindist(IntList *set1, IntList *set2)
	cdef double *obsstat_bonded(Observable_stat *stat, int j)

cdef extern from "../src/energy.h":
	cdef Observable_stat total_energy

cdef extern from "../src/energy.h":
	cdef void master_energy_calc()
	cdef void init_energies(Observable_stat *stat)


