# For C-extern Analysis 

from espresso cimport *
cimport numpy as np
from utils cimport *

cdef extern from "../src/statistics.h":
	cdef double mindist(IntList *set1, IntList *set2)
	
cdef extern from "statistics.h":
	ctypedef struct Observable_stat:
		pass
cdef extern from "energy.h":
	extern Observable_stat energy
	extern Observable_stat total_energy
