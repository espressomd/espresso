#ifndef UWERR_H
#define UWERR_H

#include <tcl.h>

/** The C implementation of the tcl function uwerr \ref tcl_uwerr.
*/
int uwerr(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
#endif
