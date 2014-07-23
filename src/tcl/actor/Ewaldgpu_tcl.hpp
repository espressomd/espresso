#ifndef EWALDGPU_TCL_H
#define EWALDGPU_TCL_H

#include "parser.hpp"
#include "config.hpp"

#ifdef EWALD_GPU

#ifndef ELECTROSTATICS
#error EWALD_GPU requires ELECTROSTATICS
#endif

int tclprint_to_result_ewaldgpu(Tcl_Interp *interp);
int tclcommand_inter_coulomb_parse_ewaldgpu(Tcl_Interp *interp, int argc, char **argv);
int tclcommand_inter_coulomb_parse_ewaldgpu_tune(Tcl_Interp * interp, int argc, char ** argv, int adaptive);
int tclcommand_inter_coulomb_parse_ewaldgpu_tunealpha(Tcl_Interp * interp, int argc, char ** argv);
#endif
#endif

