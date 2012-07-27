# For C-extern Analysis 

from espresso cimport *
cimport numpy as np
from utils cimport *

cdef extern from "../src/utils.h":
	ctypedef struct IntList:
		int *e
		int n
		int max


