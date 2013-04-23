

#ifdef __cplusplus
extern "C" {
#endif
//#include <cuda.h>
//#include <cuda_runtime.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include <stdio.h>
#include "config.h" //this is required so that the ifdefs are actually defined


#ifndef CUDA_COMMON_H
#define CUDA_COMMON_H

/** Action number for \ref mpi_get_particles. */
#define REQ_GETPARTS  16

/**cuda streams for parallel computing on cpu and gpu */
extern cudaStream_t stream[1];

extern cudaError_t err;
extern cudaError_t _err;

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

typedef struct {

  unsigned int seed;

} LB_particle_seed_gpu;

extern LB_particle_gpu *host_data;


/** This structure contains global variables associated with all of the particles and not with one individual particle */
typedef struct {
  
  /**  This is for seeding the particles' individual seeds and is initialized using irandom, beware if using for other purposes */
  unsigned int seed;
  
  unsigned int number_of_particles; 
  
  /** a boolean variable to indicate if particle info should be communicated between the cpu and gpu */
  unsigned int communication_enabled;
} GPU_global_part_vars;


void copy_forces_from_GPU();
GPU_global_part_vars* gpu_get_global_particle_vars_pointer();
LB_particle_gpu* gpu_get_particle_pointer();
LB_particle_force_gpu* gpu_get_particle_force_pointer();
LB_particle_seed_gpu* gpu_get_particle_seed_pointer();
void gpu_init_particle_comm();
void mpi_get_particles_lb(LB_particle_gpu *host_result);
void copy_part_data_to_gpu();
void mpi_send_forces_lb(LB_particle_force_gpu *host_forces);

/**cuda streams for parallel computing on cpu and gpu */
//cudaStream_t stream[1];

//cudaError_t err;
//cudaError_t _err;

/**erroroutput for memory allocation and memory copy 
 * @param err cuda error code
 * @param *file .cu file were the error took place
 * @param line line of the file were the error took place
*/

static void _cuda_safe_mem(cudaError_t err, char *file, unsigned int line){
    if( cudaSuccess != err) {                                             
      fprintf(stderr, "Cuda Memory error at %s:%u.\n", file, line);
      printf("CUDA error: %s\n", cudaGetErrorString(err));
      if ( err == cudaErrorInvalidValue )
        fprintf(stderr, "You may have tried to allocate zero memory at %s:%u.\n", file, line);
      exit(EXIT_FAILURE);
    } else {
      _err=cudaGetLastError();
      if (_err != cudaSuccess) {
        fprintf(stderr, "Error found during memory operation. Possibly however from an failed operation before. %s:%u.\n", file, line);
        printf("CUDA error: %s\n", cudaGetErrorString(err));
        if ( _err == cudaErrorInvalidValue )
          fprintf(stderr, "You may have tried to allocate zero memory before %s:%u.\n", file, line);
        exit(EXIT_FAILURE);
      }
    }
}

#define cuda_safe_mem(a) _cuda_safe_mem((a), __FILE__, __LINE__)


#define KERNELCALL(_f, _a, _b, _params) \
_f<<<_a, _b, 0, stream[0]>>>_params; \
_err=cudaGetLastError(); \
if (_err!=cudaSuccess){ \
  printf("CUDA error: %s\n", cudaGetErrorString(_err)); \
  fprintf(stderr, "error calling %s with dim %d %d %d #thpb %d in %s:%u\n", #_f, _a.x, _a.y, _a.z, _b, __FILE__, __LINE__); \
  exit(EXIT_FAILURE); \
}



#endif /* ifdef CUDA_COMMON_H */

#ifdef __cplusplus
}
#endif
