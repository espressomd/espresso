#ifndef DATAFIELD_H
#define DATAFIELD_H

#include <tcl.h>

/*************************
 * field access from tcl *
 *************************/ 

/* possible field types */
#define TYPE_INT    0
#define TYPE_DOUBLE 1

/* set callback procedure */
typedef int (SetCallback)(Tcl_Interp *interp, void *data);

typedef struct {
  void       *data;
  int         type;
  int         dimension;
  const char *name;
  SetCallback *changeproc;
} Tcl_Datafield;

/* declare the commands needed to modify particle data */
void tcl_datafield_init(Tcl_Interp *interp);

#endif
