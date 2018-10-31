#ifndef LBGPU_CUH
#define LBGPU_CUH

#include "curand_wrapper.hpp"

/** Data structure holding the velocity densities for the Lattice Boltzmann
 * system. */
typedef struct {

  /** velocity density of the node */
  float *vd;
  /** seed for the random gen */
  unsigned int *seed;
  /** flag indicating whether this site belongs to a boundary */
  unsigned int *boundary;
  curandStatePhilox4_32_10_t *philox_state;

} LB_nodes_gpu;

#endif
