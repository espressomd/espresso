from __future__ import print_function
import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

import espresso as es
import global_variables as g
try:
  es.glob.time_step=-0.01
except ValueError:
  print("Espresso does not like negative timesteps")



