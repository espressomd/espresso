
#ifndef EXTERNAL_POTENTIAL_H
#define EXTERNAL_POTENTIAL_H

#include "lattice.h"

void preinit_external_potentials();

typedef struct {
  int dummy;
} ExternalPotential;

extern ExternalPotential* external_potentials;
extern int n_external_potentials;


#endif

