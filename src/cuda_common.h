#ifdef __cplusplus
extern "C" {
#endif

#include "config.h"


#ifndef CUDA_COMMON_H
#define CUDA_COMMON_H

typedef struct {
  /** force on the particle given to md part */
  float f[3];

} LB_particle_force_gpu;

typedef struct {
  /** particle position given from md part*/
  float p[3];
  /** particle momentum struct velocity p.m->v*/
  float v[3];
#ifdef LB_ELECTROHYDRODYNAMICS
  float mu_E[3];
#endif
#ifdef ELECTROSTATICS
  float q;
#endif
  unsigned int fixed;

} LB_particle_gpu;


void cuda_enable_particle_communication();
void lb_copy_forces_GPU(LB_particle_force_gpu *host_forces);
void lb_get_particle_pointer(LB_particle_gpu** pointeradress);
void lb_get_particle_force_pointer(LB_particle_force_gpu** pointeradress);

#endif /* ifdef CUDA_COMMON_H */

#ifdef __cplusplus
}
#endif
