#ifndef IMD_H__
#define IMD_H__
/** \file imd.h 
    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/

#include <tcl.h>

int imd(ClientData data, Tcl_Interp *interp,
	int argc, char **argv);

extern int transfer_rate;

#endif
