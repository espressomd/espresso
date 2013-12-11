import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

import espresso as es
import numpy
import code_info

es.part[0].pos = (0.,0.,0.) 
es.part[0].q = 5.0
print es.part[0].q

