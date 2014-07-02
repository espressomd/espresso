#include "tcl.h"
#include "core/actor/DipolarDirectSum.hpp"
#include "actor/DipolarDirectSum.hpp"


int tclcommand_inter_magnetic_parse_dds_gpu(Tcl_Interp * interp, int argc, char ** argv)
{
 activate_dipolar_direct_sum_gpu();
  return TCL_OK;
}

int  tclprint_to_result_dds_gpu(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, " dds-gpu ", (char *) NULL);
  return TCL_OK;
}
