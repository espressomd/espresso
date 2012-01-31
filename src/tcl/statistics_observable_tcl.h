
#ifndef STATISTICS_OBSERVABLE_TCL_H
#define STATISTICS_OBSERVABLE_TCL_H

#include "config.h"
#include <tcl.h>
#include "statistics_observable.h"

int tclcommand_observable(ClientData data, Tcl_Interp *interp, int argc, char **argv);
int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ); 

//int parse_observable(Tcl_Interp* interp, int argc, char** argv, int* change, int (**A_fun)  ( void* A_args, double* A, unsigned int dim_A), int* dim_A, void** A_args);

#endif
