
#ifndef EXTERNAL_POTENTIAL_H
#define EXTERNAL_POTENTIAL_H

#include "lattice.hpp"

#define MAX_FILENAME_SIZE 32

#define EXTERNAL_POTENTIAL_TYPE_TABULATED 0
#define EXTERNAL_POTENTIAL_TYPE_ROD 1


void external_potential_pre_init();

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

void external_potential_init_energies();

int generate_external_potential(ExternalPotential** externalPotential);

int external_potential_init(int number, char* filename, int n_particle_types, double* scale);

int external_potential_tabulated_read_potential_file(int number);

void add_external_potential_forces(Particle* p);
void add_external_potential_energy(Particle* p);

int write_local_lattice_to_file(const char* filename_prefix, Lattice* lattice); 


#endif

