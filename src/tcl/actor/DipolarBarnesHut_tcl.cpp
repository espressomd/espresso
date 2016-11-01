#include "tcl.h"
#include "config.hpp"
#include "actor/DipolarBarnesHut.hpp"


#ifdef BARNES_HUT

int tclcommand_inter_magnetic_parse_gpu_bh(Tcl_Interp * interp, int argc, char ** argv)
{
 activate_dipolar_barnes_hut();
  return TCL_OK;
}

#endif

