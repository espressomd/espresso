

#include "external_potential.h"

ExternalPotential* external_potentials;
int n_external_potentials;

void preinit_external_potentials() {
  external_potentials = NULL;
  n_external_potentials = 0;
}


ExternalPotential* generate_external_potential() {
  realloc(external_potentials, (n_external_potentials+1)*sizeof(ExternalPotential));
  n_external_potentials++;

  // Test code

  Lattice l;
  double agrid[3];
  agrid[0] = agrid[1] = agrid[2] = 1;
  init_lattice(&l, agrid, 0, 0, 0);

  double* content;
  int index[3];

  lattice_allocate_memory(&l, sizeof(double));

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        double d=i*j*k;
        index[0] = i;
        index[1] = j;
        index[2] = k;
        lattice_set_data_for_local_grid_index(&l, index, &d);
      }
    }
  }
  
  for (int i=0; i<5; i++) {
    for (int j=0; j<5; j++) { 
      for (int k=0; k<5; k++) {
        double d;
        index[0] = i;
        index[1] = j;
        index[2] = k;
        lattice_get_data_for_local_halo_grid_index(&l, index, &d);
        printf("%d %d %d --> %f (%f)\n", i, j, k, d, (i-1)*(j-1)*(k-1)); 
      }
    }
  }

  

  return &external_potentials[n_external_potentials-1];
}


