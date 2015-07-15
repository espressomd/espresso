
#ifndef _IBM_TRIEL_TCL_H
#define _IBM_TRIEL_TCL_H

#include "interaction_data.hpp"

#ifdef IMMERSED_BOUNDARY

int tclcommand_inter_parse_ibm_triel(Tcl_Interp *interp, int bond_type, int argc, char **argv);
int tclprint_to_result_ibm_triel(Tcl_Interp *interp, Bonded_ia_parameters *params);

#endif

#endif
