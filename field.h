#ifndef FIELD_H
#define FIELD_H

#include <tcl.h>

/**********************************
 * standardized field declaration *
 **********************************/ 

/* possible field types */
#define TYPE_INT    0
#define TYPE_DOUBLE 1

/* possible field dimensions */
#define DIM_NPARTICLES 0
#define DIM_NCHAINS    1
#define DIM_SCALAR     2

/* readonly / read write */
#define RONLY 0
#define WRITE 1

/* callback procedure on variable set
 * the argument is the new value, the
 * value is NOT changed
 */
typedef void (VarChangeProc)(Tcl_Interp *interp);
/* actual call interface for int */
typedef void (IntVarChangeProc)(Tcl_Interp *interp, int new_val);
/* actual call interface for double */
typedef void (DoubleVarChangeProc)(Tcl_Interp *interp, double new_val);
#define NO_CALLBACK ((VarChangeProc *)NULL)

typedef struct {
  void *data;
  int  type;
  int  dimension;
  const char *name;
  VarChangeProc *changeproc;
} Field;

/* declare the commands needed to modify particle data */
void field_init(Tcl_Interp *interp);
/* init variables */
void field_var_init();
/* realloc npart */
void field_var_resize(int new_size, int dim);

#endif
