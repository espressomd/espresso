#include <tcl.h>
#ifdef USE_TK
#include <tk.h>
#endif
#include "initialize.h"

int appinit(Tcl_Interp *interp)
{
  if (Tcl_Init(interp) == TCL_ERROR)
    return (TCL_ERROR);
#ifdef USE_TK
  if (Tk_Init(interp) == TCL_ERROR)
    return (TCL_ERROR);
#endif
  if (initialize(interp) == TCL_ERROR)
    return (TCL_ERROR);
  return (TCL_OK);
}

int main(int argc, char **argv)
{
#ifdef USE_TK
  Tk_Main(argc, argv, appinit);
#else
  Tcl_Main(argc, argv, appinit);
#endif
  return 0;
}
