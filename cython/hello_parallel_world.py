
import espresso
es = espresso.espresso_instance()

if es.this_node == 0:
  N=100
### This is ugly
  tclscript="""
  setmd time_step 0.01
  setmd box_l 10. 10. 10.
  setmd skin 0.5
  thermostat langevin 1. 1.
  for { set i 0 } { $i < 100 } { incr i } {
    part $i pos [ expr 10* [ t_random ] ]  [ expr 10* [ t_random ] ]  [ expr 10* [ t_random ] ] 
  }
  integrate 1000
  """
### end ugly
  es.Tcl_Eval(tclscript)
  for i in range(N):
    es.Tcl_Eval("integrate 100")
    print "kinetic energy " + es.Tcl_Eval("analyze energy kinetic")
  es.die()
