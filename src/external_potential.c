

#include "external_potential.h"
#include "lattice.h"

ExternalPotential* external_potentials;
int n_external_potentials;

void preinit_external_potentials() {
  external_potentials = NULL;
  n_external_potentials = 0;
}


int generate_external_potential(ExternalPotential** e) {
  realloc(external_potentials, (n_external_potentials+1)*sizeof(ExternalPotential));
  n_external_potentials++;

  // Test code

  Lattice l;
  double agrid[3];
  agrid[0] = agrid[1] = agrid[2] = 1;
  int error = init_lattice(&l, agrid, 1, 0, sizeof(double));
  if (error == ES_ERROR)
      return ES_ERROR;

  double* content;
  index_t index[3];

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        double d=(i+1)*(j+1)*(k+1);
        index[0] = i;
        index[1] = j;
        index[2] = k;
        lattice_set_data_for_local_grid_index(&l, index, &d);
      }
    }
  }
  
  for (int i=1; i<4; i++) {
    for (int j=1; j<4; j++) { 
      for (int k=1; k<4; k++) {
        double *d;
        index[0] = i;
        index[1] = j;
        index[2] = k;
        lattice_get_data_for_local_halo_grid_index(&l, index, (void**) &d);
        printf("%d %d %d --> %lf (%d)\n", i, j, k, *d, (i)*(j)*(k)); 
      }
    }
  }

  

//  e = &external_potentials[n_external_potentials-1];
  return ES_OK;
}


