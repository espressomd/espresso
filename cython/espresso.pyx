


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

cdef mpi_init_helper():
  cdef int i
  cdef char** c
  i=0
  c=NULL
  mpi_init(&i, &c)

cdef extern from "../src/initialize.h":
  void on_program_start(Tcl_Interp*)


## Now comes the real deal
cdef class espresso_instance:
  cdef Tcl_Interp* interp
  def __init__(self):
    global instance_counter
    if instance_counter >= 1:
      raise Exception("Espresso shall only be instanciated once!")
    else:
      instance_counter+=1
      mpi_init_helper()
      if this_node==0:
        self.interp = Tcl_CreateInterp() 
        self.Tcl_Eval('global argv; set argv ""')
        self.Tcl_Eval('set tcl_interactive 0')
        on_program_start(self.interp)
      else:
        on_program_start(NULL)

  def __del__(self):
    raise Exception("Espresso can not be deleted")

  def Tcl_Eval(self, string):
    result=Tcl_Eval(self.interp, string)
    if result:
      raise Exception("Tcl reports an error", self.interp.result)
    return self.interp.result


## Here we create something to handle particles
cdef extern from "../src/particle_data.h":
#  int place_particle(int part, double p[3])
  pass
  

cdef class ParticleHandle:
  cdef readonly int id
  cdef bint valid
#  property position:
#    def __set__(double x[3]):
#      place_particle(id, x)
#
#      
#  def __init__(self, id):
#    self.id=id
#   

