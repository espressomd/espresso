import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

import espresso as es
import numpy
import code_info

print code_info.electrostatics_defined()
exit()

es._espressoHandle.Tcl_Eval("thermostat langevin 1. 1.")

N=100
es.glob.time_step=0.01
es.glob.skin=1.
es.glob.box_l=[10., 10., 10.]

for i in range(N):
  es.part[i].pos=numpy.random.random(3)*es.glob.box_l

es.inter[0,0].lennardJones = {"eps":1,"sigma":1,"shift":0.25}

es._espressoHandle.Tcl_Eval("integrate 100")
for i in range(N):
  print es.part[i].pos

es._espressoHandle.die()

