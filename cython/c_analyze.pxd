# For C-extern Analysis 

from espresso cimport *
cimport numpy as np
from utils cimport *

cdef extern from "../src/statistics.h":
	cdef double mindist(IntList *set1, IntList *set2)
	
cdef extern from "../src/statistics.h":
	ctypedef struct Observable_stat:
		int init_status
		DoubleList data

cdef extern from "../src/energy.h":
	cdef Observable_stat energy
	cdef Observable_stat total_energy

cdef extern from "../src/energy.h":
	cdef void master_energy_calc()
	cdef void init_energies(Observable_stat *stat)
