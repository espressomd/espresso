#include "global.h"
#include "field.h"

/******************* variables ****************/

int N;
double  *x,  *y,  *z;
double  *q;
int     *ident;
double *vx, *vy, *vz;
double *Fx, *Fy, *Fz;

/*********** field interface data *************/
Field fields[] = {
  { (void *)&N, TYPE_INT, DIM_SCALAR, "N",
    (VarChangeProc *)set_npart },
  { (void *)&Nchains, TYPE_INT, DIM_SCALAR, "Nchains",
    (VarChangeProc *)set_nchains },
  { (void *)&x, TYPE_DOUBLE, DIM_NPARTICLES, "x", NO_CALLBACK },
  { (void *)&y, TYPE_DOUBLE, DIM_NPARTICLES, "y", NO_CALLBACK },
  { (void *)&z, TYPE_DOUBLE, DIM_NPARTICLES, "z", NO_CALLBACK },
  { (void *)&q, TYPE_DOUBLE, DIM_NPARTICLES, "q", NO_CALLBACK },
  { (void *)&ident, TYPE_INT, DIM_NPARTICLES, "ident", NO_CALLBACK },
  { (void *)&vx, TYPE_DOUBLE, DIM_NPARTICLES, "vx", NO_CALLBACK },
  { (void *)&vy, TYPE_DOUBLE, DIM_NPARTICLES, "vy", NO_CALLBACK },
  { (void *)&vz, TYPE_DOUBLE, DIM_NPARTICLES, "vz", NO_CALLBACK },
  { (void *)&Fx, TYPE_DOUBLE, DIM_NPARTICLES, "fx", NO_CALLBACK },
  { (void *)&Fy, TYPE_DOUBLE, DIM_NPARTICLES, "fy", NO_CALLBACK },
  { (void *)&Fz, TYPE_DOUBLE, DIM_NPARTICLES, "fz", NO_CALLBACK },
  { NULL, -1, -1, "DO NOT DELETE THIS ENTRY!", NO_CALLBACK }
};

/***************** procedures ******************/

void init_variables()
{
  field_var_init();
}

void set_npart(Tcl_Interp *interp, int npart)
{
  field_var_resize(npart, DIM_NPARTICLES);

  N = npart;
}

void set_nchains(Tcl_Interp *interp, int nchains)
{
  field_var_resize(nchains, DIM_NCHAINS);
  chain_set_nchains(nchains);

  Nchains = nchains;
}
