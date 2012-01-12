import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

import espresso as es
import numpy

es._espressoHandle.Tcl_Eval("thermostat lb 1.")

N=100
es.glob.time_step=0.01
es.glob.skin=1.
es.glob.box_l=[10., 10., 10.]
es.lb["cpu"].agrid=1
es.lb["cpu"].dens=1
es.lb["cpu"].friction=1
es.lb["cpu"].tau=0.1
es.lb["cpu"].visc=1

print es.lb["cpu"].dens
print es.lb["cpu"].visc
print es.lb["cpu"].agrid
es.lb["cpu"].print_vtk_velocity="test.vtk"
for i in range(N):
  es.part[i].pos=numpy.random.random(3)*es.glob.box_l

es.inter[0,0].lennardJones = {"eps":1,"sigma":1,"shift":0.25}

es._espressoHandle.Tcl_Eval("integrate 100")
#for i in range(N):
#  print es.part[i].pos

es._espressoHandle.die()

