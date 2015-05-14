
#ifndef _AFFINITY_TCL_H
#define _AFFINITY_TCL_H

#include "../parser.hpp"

#ifdef AFFINITY

int tclprint_to_result_affinityIA(Tcl_Interp *interp, int i, int j);

int tclcommand_inter_parse_affinity(Tcl_Interp * interp,
				int part_type_a, int part_type_b,
				int argc, char ** argv);

#endif

#endif
