import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

import espresso as es
import numpy

es._espressoHandle.Tcl_Eval("thermostat lb 1.")
dev="cpu"
N=100
es.glob.time_step=0.01
es.glob.skin=1.
es.glob.box_l=[10., 10., 10.]
es.cu.device=0
print "cuda device:"
print es.cu.device
es.lb[dev].agrid=1
es.lb[dev].dens=1
es.lb[dev].visc=1
es.lb[dev].friction=1
es.lb[dev].tau=0.1

es.lb[dev].ext_force=[1., 2., 3.,]
print es.lb[dev].ext_force
print es.lb[dev].dens
print es.lb[dev].visc
print es.lb[dev].agrid
es.lb[dev].print_vtk_velocity="test.vtk"
#es.lb[dev].checkpoint_style=1
#es.lb[dev].checkpoint="cp.dat"
for i in range(N):
  es.part[i].pos=numpy.random.random(3)*es.glob.box_l

es.inter[0,0].lennardJones = {"eps":1,"sigma":1,"shift":0.25}

es._espressoHandle.Tcl_Eval("integrate 100")
#for i in range(N):
#  print es.part[i].pos

es._espressoHandle.die()

