#ifndef BENDING_FORCE_IBM_TCL_H
#define BENDING_FORCE_IBM_TCL_H

#include "parser.hpp"
#include "interaction_data.hpp"

#ifdef BENDING_FORCE_IMMERSED_BOUNDARY

int tclcommand_inter_parse_bending_force_ibm(Tcl_Interp *interp, int bond_type, int argc, char **argv);
int tclprint_to_result_bending_force_ibmIA(Tcl_Interp *interp, Bonded_ia_parameters *params);

#endif
#endif
