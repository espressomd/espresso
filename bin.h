#ifndef BIN_H
#define BIN_H

#include "utils.h"
#include <tcl.h>

/** Implementation of the tcl command bin, which can be used
    to bin data into arbitrary bins
*/
int bin(ClientData data, Tcl_Interp *interp,
	int argc, char **argv);

#endif
