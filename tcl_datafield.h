#ifndef DATAFIELD_H
#define DATAFIELD_H

#include <tcl.h>

/*************************
 * field access from tcl *
 *************************/ 

/* declare the commands needed to modify particle data, set
   interaction parameters and global parameters, i. e.
   part, inter and setmd */
void tcl_datafield_init(Tcl_Interp *interp);

#endif
