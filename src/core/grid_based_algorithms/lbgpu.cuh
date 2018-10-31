#ifndef LBGPU_CUH
#define LBGPU_CUH

#include "config.hpp"

#ifdef CUDA
#include "curand_wrapper.hpp"

#ifdef LB_GPU
/** Data structure holding the velocity densities for the Lattice Boltzmann
 * system. */
typedef struct {

  /** velocity density of the node */
  float *vd;
  /** seed for the random gen */
  unsigned int *seed;
  /** flag indicating whether this site belongs to a boundary */
  unsigned int *boundary;

} LB_nodes_gpu;
#endif // LB_GPU

#endif //CUDA
#endif
