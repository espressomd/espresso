#ifndef CHAIN_H
#define CHAIN_H

typedef struct {
  int start, end;
  int ident;
} Chain;

#include <tcl.h>

/* declare the commands to define the chains */
void chain_init(Tcl_Interp *interp);

/* callback on nchains set */
void chain_set_nchains(int newVal);

#endif
