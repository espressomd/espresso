#ifndef IMD_H__
#define IMD_H__

#include <tcl.h>

int imd(ClientData data, Tcl_Interp *interp,
	int argc, char **argv);

extern int transfer_rate;

#endif
