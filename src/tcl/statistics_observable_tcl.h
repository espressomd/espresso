
#ifndef STATISTICS_OBSERVABLE_TCL_H
#define STATISTICS_OBSERVABLE_TCL_H

#include "config.h"
#include <tcl.h>
#include "statistics_observable.h"

int tclcommand_observable(ClientData data, Tcl_Interp *interp, int argc, char **argv);
int tclcommand_observable_print_formatted(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs, double* values);
int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ); 

#endif
