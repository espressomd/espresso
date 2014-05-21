#ifndef MMM1D_GPU_TCL_H
#define MMM1D_GPU_TCL_H

#include "parser.hpp"
#include "config.hpp"

#ifdef MMM1D_GPU

#ifndef ELECTROSTATICS
#error MMM1D_GPU requires ELECTROSTATICS
#endif

/// print the mmm1d parameters to the interpreters result
int tclprint_to_result_MMM1DGPU(Tcl_Interp *interp);

/// parse the mmm1d parameters
int tclcommand_inter_coulomb_parse_mmm1dgpu(Tcl_Interp *interp, int argc, char **argv);

#endif
#endif
