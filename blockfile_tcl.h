/** \file blockfile_tcl.h
    Contains only the tcl interface for block coded files.

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    It is the header file for \ref blockfile_tcl.c "blockfile_tcl.c" and provides the
    tcl command \ref tcl_blockfile.
*/
#ifndef BLOCKFILE_TCL_H
#define BLOCKFILE_TCL_H
#include <tcl.h>

/** Implementation of the Tcl command \ref tcl_blockfile. Allows to read and write
    blockfile comfortably from Tcl. */
int blockfile(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv);
#endif
