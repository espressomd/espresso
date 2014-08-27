import sys

## Let's import on_program_start
cdef extern from "communication.hpp":
  void mpi_init(int* argc = NULL, char ***argv = NULL)
  int this_node
  void mpi_stop()

cdef extern from "initialize.hpp":
  void on_program_start()
  void mpi_loop()

### Here we make a minimalistic Tcl_Interp available
# cdef extern from "tcl.h":
#   cdef struct Tcl_Interp:
#     pass

#   Tcl_Interp* Tcl_CreateInterp()
#   int Tcl_Eval(Tcl_Interp * interp, char* script)
#   char* Tcl_GetStringResult(Tcl_Interp * interp)

# cdef extern from "tcl/initialize_interpreter.hpp":
#   int appinit(Tcl_Interp *interp)
  
# cdef class TclInterpreter:
#     # Define Tcl interpreter
#     cdef Tcl_Interp* interp

#     def __init__(self):
#         self.interp = Tcl_CreateInterp()
#         appinit(self.interp)
    
#     def eval(self, string):
#         result = Tcl_Eval(self.interp, string)
#         tclresult = Tcl_GetStringResult(self.interp)
#         if result:
#             raise Exception("Tcl reports an error", tclresult)
  
#         return tclresult

def setup():
    # Main code
    mpi_init()

    # Main slave loop
    if this_node != 0:
        on_program_start()
        mpi_loop()
        sys.exit()

    # # Initialize Tcl
    # tcl = TclInterpreter()
    # tcl.eval('global argv; set argv ""')
    # tcl.eval('set tcl_interactive 0')
    on_program_start()

    # return tcl

def finalize():
    mpi_stop()
