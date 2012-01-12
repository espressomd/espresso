cimport numpy as np
import numpy as np
import particle_data
cimport particle_data 
import interaction_data
cimport interaction_data
import global_variables
from integrate import integrate
import thermostat
cimport thermostat
from changeVolume import changeVolume
from invalidateSystem import invalidateSystem
import cellsystem
cimport cellsystem

import debye_hueckel
#cimport myconfig
#import utils

#public enum:
#  ERROR=-1


#### This is the main Espresso Cython interface.

### Now come all the includes

cdef int instance_counter 

### Here we make a minimalistic Tcl_Interp available
cdef extern from "tcl.h":
  cdef struct Tcl_Interp:
    char *result
    int errorLine
    
  Tcl_Interp* Tcl_CreateInterp()
  int Tcl_Eval (Tcl_Interp * interp, char* script)
  ctypedef Tcl_Interp* Interp_pointer


## Let's import on_program_start
cdef extern from *:
  void mpi_init(int* argc, char*** argv)
  int this_node
  void mpi_stop()

cdef mpi_init_helper():
  cdef int i
  cdef char** c
  i=0
  c=NULL
  mpi_init(&i, &c)

cdef extern from "initialize.h":
  void on_program_start(Tcl_Interp*)
  void mpi_loop()


## Now comes the real deal
cdef class EspressoHandle:
  cdef Tcl_Interp* interp
  cdef public int this_node
  def __init__(self):
    global instance_counter
    if instance_counter >= 1:
      raise Exception("Espresso shall only be instanciated once!")
    else:
      instance_counter+=1
      mpi_init_helper()
      self.this_node=this_node
      if this_node==0:
        self.interp = Tcl_CreateInterp() 
        self.Tcl_Eval('global argv; set argv ""')
        self.Tcl_Eval('set tcl_interactive 0')
        on_program_start(self.interp)
      else:
        on_program_start(NULL)
        mpi_loop()
        self.die()

  def __del__(self):
    self.die()
    raise Exception("Espresso can not be deleted")

  def Tcl_Eval(self, string):
    result=Tcl_Eval(self.interp, string)
    if result:
      raise Exception("Tcl reports an error", self.interp.result)
    return self.interp.result

  def die(self):
    mpi_stop()

_espressoHandle=EspressoHandle()
glob=global_variables.GlobalsHandle()
part=particle_data.particleList()

inter=interaction_data.InteractionList()

if this_node==0:
  glob=global_variables.GlobalsHandle()
else:
  exit()

      

#
#      
#  def __init__(self, id):
#    self.id=id
#   





