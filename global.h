#ifndef GLOBAL_H
#define GLOBAL_H
#include "field.h"
#include "chain.h"

/**********************************************
 * global variables
 *
 * init_variables, set_npart and set_nchains
 * can be used to keep track of the allocation
 *
 * for simple fields an automatic mechanism is
 * provided by field.c/h.
 **********************************************/


/******************* variables ****************/
/***** particle data ******/
extern int N;
/* position */
extern double  *x,  *y,  *z;
/* charge */
extern double  *q;
/* ident */
extern int *ident;
/* velocity */
extern double *vx, *vy, *vz;
/* force */
extern double *Fx, *Fy, *Fz;

/***** group data *********
 * chain.c                */
extern int Nchains;
extern Chain chains[];

/*********** field interface data ************/
extern Field fields[];

/****************** procedures ***************/
/* initialize variables */
void init_variables();
/* reallocate npart */
void set_npart(Tcl_Interp *interp, int npart);
/* reallocate nchains */
void set_nchains(Tcl_Interp *interp, int npart);

#endif
