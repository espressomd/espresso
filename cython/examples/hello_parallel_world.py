
import espresso as es
import particle_data
import numpy as np
import espresso.glob as gl

if es.this_node == 0:
  N=1
  es.time_step=0.01
  es.skin=1.
  es.box_l=[10., 10., 10.]

  es.Tcl_Eval("thermostat langevin 1. 1.")

#  es.part[0].pos=[0., 0., 0.]
#  for p in es.part:
#    print p.pos
#
##  inter[0,0]=es.inter.lennard-jones(sigma=1, epsilon=1)
##  es.part[:].pos+=[1.,0., 0.]
##

  es.Tcl_Eval("integrate 100")
  for i in range(N):
    p=particle_data.ParticleHandle(i)
    print p.pos

### This is ugly




#  tclscript="""
#  setmd time_step 0.01
#  setmd box_l 10. 10. 10.
#  setmd skin 0.5
#  thermostat langevin 1. 1.
#  for { set i 0 } { $i < 100 } { incr i } {
#    part $i pos [ expr 10* [ t_random ] ]  [ expr 10* [ t_random ] ]  [ expr 10* [ t_random ] ] 
#  }
#  integrate 1000
#  """
#### end ugly
#  es.Tcl_Eval(tclscript)
#  p=espresso.ParticleHandle(0)
#  print p.pos
#  print es.Tcl_Eval("part 0 print pos")
#  exit()
#  for i in range(N):
#    es.Tcl_Eval("integrate 100")
#    print "kinetic energy " + es.Tcl_Eval("analyze energy kinetic")
#  es.die()
