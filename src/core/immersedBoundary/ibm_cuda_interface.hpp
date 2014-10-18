

// *******
// This is an internal file of the IMMERSED BOUNDARY implementation
// It should not be included by any main Espresso routines
// Functions to be exported for Espresso are in ibm_main.hpp

#ifndef IBM_CUDA_INTERFACE_HPP
#define IBM_CUDA_INTERFACE_HPP

#ifdef IMMERSED_BOUNDARY

// *********** Communication functions ********
// Implemented in real C++, but called from the ibm_cuda.cu
void IBM_cuda_mpi_send_velocities();
void IBM_cuda_mpi_get_particles();

void ParticleVelocitiesFromLB_GPU();

// ******** data types for CUDA and MPI communication ******
typedef struct {
  float pos[3];
  float f[3];
  bool isVirtual;
} IBM_CUDA_ParticleDataInput;

typedef struct {
  float v[3];
} IBM_CUDA_ParticleDataOutput;

// ******** global variables for CUDA and MPI communication ******
extern IBM_CUDA_ParticleDataInput *IBM_ParticleDataInput_host;
extern IBM_CUDA_ParticleDataOutput *IBM_ParticleDataOutput_host;



#endif

#endif