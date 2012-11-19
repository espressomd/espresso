
#ifndef EXTERNAL_POTENTIAL_H
#define EXTERNAL_POTENTIAL_H

#include "lattice.h"

#define MAX_FILENAME_SIZE 32

#define EXTERNAL_POTENTIAL_TYPE_TABULATED 0
#define EXTERNAL_POTENTIAL_TYPE_ROD 1


void preinit_external_potentials();

typedef struct {
  char filename[MAX_FILENAME_SIZE];
  Lattice potential;
} ExternalPotentialTabulated;

typedef struct {
  int dummy;
  int type;
  double* scale;
  int n_particle_types;
  union {
    ExternalPotentialTabulated tabulated;
  } e;
  double energy;
} ExternalPotential;

extern ExternalPotential* external_potentials;
extern int n_external_potentials;

int external_potential_tabulated_init(int number, char* filename, int n_particle_types, double* scale);

int generate_external_potential(ExternalPotential** externalPotential);

int external_potential_init(int number, char* filename, int n_particle_types, double* scale);

int external_potential_tabulated_read_potential_file(int number);

void add_external_potential_forces(Particle* p);

#endif

