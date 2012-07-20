# For C-extern Analysis 

from espresso cimport *
cimport numpy as np
from utils cimport *

cdef extern from "../src/statistics.h":
	cdef double mindist(IntList *set1, IntList *set2)
