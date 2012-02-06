#ifndef INITIALIZE_INTERPRETER_H
#define INITIALIZE_INTERPRETER_H
#include "config.h"
#include <tcl.h>

void register_tcl_commands(Tcl_Interp* interp);

#endif
