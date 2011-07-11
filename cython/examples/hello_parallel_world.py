
import espresso as es
import numpy

es._espressoHandle.Tcl_Eval("thermostat langevin 1. 1.")

N=100
es.glob.time_step=0.01
es.glob.skin=1.
es.glob.box_l=[10., 10., 10.]

for i in range(N):
  p=es.particle_data.ParticleHandle(i)
  p.pos=numpy.random.random(3)*es.glob.box_l


es._espressoHandle.Tcl_Eval("integrate 100")
for i in range(N):
  p=es.particle_data.ParticleHandle(i)
  print p.pos

es._espressoHandle.die()

