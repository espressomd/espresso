#ifndef TRIEL_TCL_H
#define TRIEL_TCL_H
//Interface for triel.h/c

#include "parser.hpp"
#include "interaction_data.hpp"

#ifdef TRIELASTIC


int tclcommand_inter_parse_triel(Tcl_Interp *interp, int bond_type, int argc, char **argv);


int tclprint_to_result_trielIA(Tcl_Interp *interp, Bonded_ia_parameters *params);

#endif

#endif
