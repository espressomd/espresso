#include "actor/EwaldgpuForce.hpp"
#include "cuda_utils.hpp"

#ifdef EWALD_GPU

#include <stdio.h>
#include <iostream>

#include "interaction_data.hpp"
#include "EspressoSystemInterface.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

typedef ewaldgpu_real real;
Ewaldgpu_params ewaldgpu_params;

//Stream
cudaEvent_t *start, *stop;
cudaStream_t    *stream0;

void cuda_check_error(const dim3 &block, const dim3 &grid,const char *function, const char *file, unsigned int line);

//Error handler
static void HandleError(cudaError_t err,const char *file,int line)
{
  if( cudaSuccess != err) {
    fprintf(stderr, "Cuda Memory error at %s:%u.\n", file, line);
    printf("CUDA error: %s\n", cudaGetErrorString(err));
    if ( err == cudaErrorInvalidValue )
      fprintf(stderr, "You may have tried to allocate zero memory at %s:%u.\n", file, line);
    exit(EXIT_FAILURE);
  } else {
    err=cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "Error found during memory operation. Possibly however from an failed operation before. %s:%u.\n", file, line);
      printf("CUDA error: %s\n", cudaGetErrorString(err));
      if ( err == cudaErrorInvalidValue )
	fprintf(stderr, "You may have tried to allocate zero memory before %s:%u.\n", file, line);
      exit(EXIT_FAILURE);
    }
  }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define HANDLE_NULL( a ) {if (a == NULL) {printf( "Host memory failed in %s at line %d\n", __FILE__, __LINE__ ); exit( EXIT_FAILURE );}}

//Kernels
/*Here MaxThreads/LowThreads means that, if for example 10000 Particles/k-Vectors are given and the GPU uses 256 Threads the first 19*(2*256)=9728 will be
  reduced by the MaxThread function and the remaining 272 Particles will be reduced by the LowThread function*/
//Structure factor
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_Rho_hat_MaxThreads(real *k,real *r_i, real *q_i, real *rho_hat, int N, int num_k, int maxThreadsPerBlock,int loops)
{
  //Reject unneeded blocks
  if (blockIdx.x*gridDim.x+blockIdx.y >= num_k) return;

  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  real kr;
  real factor;
  real sin_kr;
  real cos_kr;
  real *sin_ptr=&sin_kr;
  real *cos_ptr=&cos_kr;
  int blockId = blockIdx.x*gridDim.x+blockIdx.y;
  int l = loops;

  //Init
  i = tid;
  sdata[tid] = 0;
  sdata[tid+2*blockSize] = 0;

  //Reduction
  while (i < mTPB)
    {
      kr = k[blockId]*r_i[3*(i+2*l*mTPB)]+k[blockId+num_k]*r_i[3*(i+2*l*mTPB)+1]+k[blockId+2*num_k]*r_i[3*(i+2*l*mTPB)+2];
      sincos(kr,sin_ptr,cos_ptr);
      factor = q_i[i+2*l*mTPB];
      sdata[tid]      						+=  factor*cos_kr;
      sdata[tid+2*blockSize]      += -factor*sin_kr;

      //BECAUSE nIsPow2=True
      kr = k[blockId]*r_i[3*(i+2*l*mTPB+blockSize)]+k[blockId+num_k]*r_i[3*(i+2*l*mTPB+blockSize)+1]+k[blockId+2*num_k]*r_i[3*(i+2*l*mTPB+blockSize)+2];
      sincos(kr,sin_ptr,cos_ptr);
      factor = q_i[i+2*l*mTPB+blockSize];
      sdata[tid]      						+=  factor*cos_kr;
      sdata[tid+2*blockSize]      += -factor*sin_kr;

      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 64];  }  __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      rho_hat[blockId]       += sdata[0];
      rho_hat[blockId+num_k] += sdata[2*blockSize];
    }
  __syncthreads();
}
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_Rho_hat_LowThreads(real *k,real *r_i, real *q_i, real *rho_hat, int N, int num_k,int maxThreadsPerBlock,int elapsedLoops)
{
  //Reject unneeded blocks
  if (blockIdx.x*gridDim.x+blockIdx.y >= num_k) return;

  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  int l=elapsedLoops;
  real kr;
  real factor;
  real sin_kr;
  real cos_kr;
  real *sin_ptr=&sin_kr;
  real *cos_ptr=&cos_kr;
  int blockId = blockIdx.x*gridDim.x+blockIdx.y;

  //Init
  i = tid;
  sdata[tid] = 0;
  sdata[tid+2*blockSize] = 0;

  //Reduction
  while (i < N-2*l*mTPB)
    {
      kr = k[blockId]*r_i[3*(i+2*l*mTPB)]+k[blockId+num_k]*r_i[3*(i+2*l*mTPB)+1]+k[blockId+2*num_k]*r_i[3*(i+2*l*mTPB)+2];
      sincos(kr,sin_ptr,cos_ptr);
      factor = q_i[i+2*l*mTPB];

      sdata[tid]      						+=  factor*cos_kr;
      sdata[tid+2*blockSize]      += -factor*sin_kr;
      if (nIsPow2 || i + blockSize < N-2*l*mTPB)
	{
	  kr = k[blockId]*r_i[3*(i+blockSize+2*l*mTPB)]+k[blockId+num_k]*r_i[3*(i+blockSize+2*l*mTPB)+1]+k[blockId+2*num_k]*r_i[3*(i+blockSize+2*l*mTPB)+2];
	  sincos(kr,sin_ptr,cos_ptr);
	  factor = q_i[i+2*l*mTPB+blockSize];
	  sdata[tid]      						+=  factor*cos_kr;
	  sdata[tid+2*blockSize]      += -factor*sin_kr;
	}
      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 64];  }  __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      rho_hat[blockId]       += sdata[0];
      rho_hat[blockId+num_k] += sdata[2*blockSize];
    }
}
//Forces in reciprocal space
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_ForcesReci_MaxThreads(real *k,real *r_i, real *q_i, real *rho_hat, real *infl_factor, real *forces_reci, int N, int num_k, real V, real coulomb_prefactor, int maxThreadsPerBlock,int loops)
{
  //Reject unneeded blocks
  if (blockIdx.x*gridDim.x+blockIdx.y >= N) return;

  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  real kr;
  real factor;
  real sin_kr;
  real cos_kr;
  real *sin_ptr=&sin_kr;
  real *cos_ptr=&cos_kr;
  int blockId = blockIdx.x*gridDim.x+blockIdx.y;
  int l = loops;

  //Init
  i = tid;
  sdata[tid] 						 = 0;
  sdata[tid+2*blockSize] = 0;
  sdata[tid+4*blockSize] = 0;

  //Reduction
  while (i < mTPB)
    {
      kr = k[i+2*l*mTPB]*r_i[3*blockId] + k[i+2*l*mTPB+num_k]*r_i[3*blockId+1] + k[i+2*l*mTPB+2*num_k]*r_i[3*blockId+2];
      sincos(kr,sin_ptr,cos_ptr);
      factor = infl_factor[i+2*l*mTPB] * q_i[blockId];

      sdata[tid]      			 += factor * (k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB]         * sin_kr + k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB+num_k]         * cos_kr);
      sdata[tid+2*blockSize] += factor * (k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB]   * sin_kr + k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB+num_k]   * cos_kr);
      sdata[tid+4*blockSize] += factor * (k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB] * sin_kr + k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB+num_k] * cos_kr);
      //BECAUSE nIsPow2=True
      kr = k[i+2*l*mTPB+blockSize]*r_i[3*blockId] + k[i+2*l*mTPB+num_k+blockSize]*r_i[3*blockId+1] + k[i+2*l*mTPB+2*num_k+blockSize]*r_i[3*blockId+2];
      factor = infl_factor[i+2*l*mTPB+blockSize]* q_i[blockId];
      sincos(kr,sin_ptr,cos_ptr);
      sdata[tid]      			 += factor * (k[i+2*l*mTPB+blockSize]*rho_hat[i+2*l*mTPB+blockSize]         * sin_kr + k[i+2*l*mTPB+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize]         * cos_kr);
      sdata[tid+2*blockSize] += factor * (k[i+2*l*mTPB+num_k+blockSize]*rho_hat[i+2*l*mTPB+blockSize]   * sin_kr + k[i+2*l*mTPB+num_k+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize]   * cos_kr);
      sdata[tid+4*blockSize] += factor * (k[i+2*l*mTPB+2*num_k+blockSize]*rho_hat[i+2*l*mTPB+blockSize] * sin_kr + k[i+2*l*mTPB+2*num_k+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize] * cos_kr);

      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 512];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 256];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 128];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 64]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 64];  } __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 32];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 16];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 8]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 4]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 2]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 1]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      forces_reci[3*blockId]     += sdata[0]*coulomb_prefactor;
      forces_reci[3*blockId+1]   += sdata[2*blockSize]*coulomb_prefactor;
      forces_reci[3*blockId+2] += sdata[4*blockSize]*coulomb_prefactor;
    }
}

template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_ForcesReci_LowThreads(real *k,real *r_i, real *q_i, real *rho_hat, real *infl_factor, real *forces_reci, int N, int num_k, real V, real coulomb_prefactor, int maxThreadsPerBlock,int elapsedLoops)
{
  //Reject unneeded blocks
  if (blockIdx.x*gridDim.x+blockIdx.y >= N) return;

  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  int l=elapsedLoops;
  real kr;
  real factor;
  real sin_kr;
  real cos_kr;
  real *sin_ptr=&sin_kr;
  real *cos_ptr=&cos_kr;
  int blockId = blockIdx.x*gridDim.x+blockIdx.y;

  //Init
  i = tid;
  sdata[tid] 						 = 0;
  sdata[tid+2*blockSize] = 0;
  sdata[tid+4*blockSize] = 0;
  while (i < num_k-2*l*mTPB)
    {
      kr = k[i+2*l*mTPB]*r_i[3*blockId] + k[i+2*l*mTPB+num_k]*r_i[3*blockId+1] + k[i+2*l*mTPB+2*num_k]*r_i[3*blockId+2];
      sincos(kr,sin_ptr,cos_ptr);
      factor = infl_factor[i+2*l*mTPB] * q_i[blockId];

      //Reduction
      sdata[tid]      			 += factor * (k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB]         * sin_kr + k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB+num_k]         * cos_kr);
      sdata[tid+2*blockSize] += factor * (k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB]   * sin_kr + k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB+num_k]   * cos_kr);
      sdata[tid+4*blockSize] += factor * (k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB] * sin_kr + k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB+num_k] * cos_kr);
      if (nIsPow2 || i + blockSize < num_k-2*l*mTPB)
	{
	  kr = k[i+2*l*mTPB+blockSize]*r_i[3*blockId] + k[i+2*l*mTPB+num_k+blockSize]*r_i[3*blockId+1] + k[i+2*l*mTPB+2*num_k+blockSize]*r_i[3*blockId+2];
	  sincos(kr,sin_ptr,cos_ptr);
	  factor = infl_factor[i+2*l*mTPB+blockSize] * q_i[blockId];
	  sdata[tid]      			 += factor * (k[i+2*l*mTPB+blockSize]*rho_hat[i+2*l*mTPB+blockSize]         * sin_kr + k[i+2*l*mTPB+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize]         * cos_kr);
	  sdata[tid+2*blockSize] += factor * (k[i+2*l*mTPB+num_k+blockSize]*rho_hat[i+2*l*mTPB+blockSize]   * sin_kr + k[i+2*l*mTPB+num_k+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize]   * cos_kr);
	  sdata[tid+4*blockSize] += factor * (k[i+2*l*mTPB+2*num_k+blockSize]*rho_hat[i+2*l*mTPB+blockSize] * sin_kr + k[i+2*l*mTPB+2*num_k+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize] * cos_kr);
	}
      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 512];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 256];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 128];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 64]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 64];  } __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 32];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 16];sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 8]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 4]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 2]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; sdata[tid+2*blockSize] += sdata[tid+2*blockSize + 1]; sdata[tid+4*blockSize] += sdata[tid+4*blockSize + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      forces_reci[3*blockId]   += sdata[0]*coulomb_prefactor;
      forces_reci[3*blockId+1] += sdata[2*blockSize]*coulomb_prefactor;
      forces_reci[3*blockId+2] += sdata[4*blockSize]*coulomb_prefactor;
    }
}

//Energy in reciprocal space
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_EnergyReci_MaxThreads(real *rho_hat, real *infl_factor,real *energy_reci, int N, int num_k, real V,int maxThreadsPerBlock,int loops,real coulomb_prefactor)
{
  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  real factor;
  int l = loops;

  //Init
  i = tid;
  sdata[tid] = 0;

  //Reduction
  while (i < mTPB)
    {
      factor = infl_factor[i+2*l*mTPB] / 2;

      sdata[tid]      			 += factor * (rho_hat[i+2*l*mTPB]*rho_hat[i+2*l*mTPB] + rho_hat[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB+num_k]);
      //BECAUSE nIsPow2=True
      factor = infl_factor[i+2*l*mTPB+blockSize] / 2;
      sdata[tid]      			 += factor * (rho_hat[i+2*l*mTPB+blockSize]*rho_hat[i+2*l*mTPB+blockSize] + rho_hat[i+2*l*mTPB+num_k+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize]);

      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64];  } __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      energy_reci[0] += sdata[0]*coulomb_prefactor;
    }
  __syncthreads();
}
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_EnergyReci_LowThreads(real *rho_hat, real *infl_factor,real *energy_reci, int N, int num_k, real V,int maxThreadsPerBlock,int elapsedLoops,real coulomb_prefactor)
{
  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  int l=elapsedLoops;
  real factor;

  //Init
  i = tid;
  sdata[tid] = 0;

  //Reduction
  while (i < num_k-2*l*mTPB)
    {
      factor = infl_factor[i+2*l*mTPB] / 2;

      sdata[tid]    += factor * (rho_hat[i+2*l*mTPB]*rho_hat[i+2*l*mTPB] + rho_hat[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB+num_k]);
      if (nIsPow2 || i + blockSize < num_k-2*l*mTPB)
	{
	  factor = infl_factor[i+2*l*mTPB+blockSize] / 2;
	  sdata[tid]  += factor * (rho_hat[i+2*l*mTPB+blockSize]*rho_hat[i+2*l*mTPB+blockSize] + rho_hat[i+2*l*mTPB+num_k+blockSize]*rho_hat[i+2*l*mTPB+num_k+blockSize]);
	}
      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64];  } __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      energy_reci[0] += sdata[0]*coulomb_prefactor;
    }
}
//q squared
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_q_sqr_MaxThreads(real *q_i, real *q_sqr, int N, int maxThreadsPerBlock,int loops)
{
  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  int l = loops;

  //Init
  i = tid;
  sdata[tid] 						 = 0;

  //Reduction
  while (i < mTPB)
    {
      sdata[tid]      			 += q_i[i+2*l*mTPB]*q_i[i+2*l*mTPB];
      //BECAUSE nIsPow2=True
      sdata[tid]      			 += q_i[i+2*l*mTPB+blockSize]*q_i[i+2*l*mTPB+blockSize];

      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64];  } __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      q_sqr[0] += sdata[0];
    }
  __syncthreads();
}
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_q_sqr_LowThreads(real *q_i, real *q_sqr, int N, int maxThreadsPerBlock,int elapsedLoops)
{
  //Variables
  extern __shared__ real sdata[];
  int tid = threadIdx.x;
  int i = tid;
  int gridSize = blockSize*2*gridDim.x;
  int mTPB=maxThreadsPerBlock;
  int l=elapsedLoops;

  //Init
  i = tid;
  sdata[tid] 						 = 0;

  //Reduction
  while (i < N-2*l*mTPB)
    {
      sdata[tid]    += q_i[i+2*l*mTPB]*q_i[i+2*l*mTPB];
      if (nIsPow2 || i + blockSize < N-2*l*mTPB)
	{
	  sdata[tid]  += q_i[i+2*l*mTPB+blockSize]*q_i[i+2*l*mTPB+blockSize];
	}
      i += gridSize;
    }
  __syncthreads();

  if (blockSize >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid + 64];  } __syncthreads(); }
  if (tid < 32) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; __syncthreads(); }
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; __syncthreads(); }
    if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; __syncthreads(); }
    if (blockSize >= 8)  {sdata[tid] += sdata[tid + 4]; __syncthreads(); }
    if (blockSize >= 4)  {sdata[tid] += sdata[tid + 2]; __syncthreads(); }
    if (blockSize >= 2)  {sdata[tid] += sdata[tid + 1]; __syncthreads(); }
  }
  //Write result for this block to global memory
  if (tid == 0)
    {
      q_sqr[0] += sdata[0];
    }
}

//Ewaldgpuforce
EwaldgpuForce::EwaldgpuForce(SystemInterface &s, double rcut, int num_kx, int num_ky, int num_kz, double alpha)
  :m_dev_k(NULL), m_k(NULL), m_dev_rho_hat(NULL), m_rho_hat(NULL), m_dev_infl_factor(NULL), m_infl_factor(NULL), m_dev_energy_reci(NULL), m_energy_reci(NULL), m_dev_q_sqr(NULL), m_q_sqr(NULL)
{
  //Interface sanity checks
  if(!s.requestFGpu())
    std::cerr << "EwaldgpuForce needs access to forces on GPU!" << std::endl;

  if(!s.requestRGpu())
    std::cerr << "EwaldgpuForce needs access to positions on GPU!" << std::endl;

  if(!s.requestQGpu())
    std::cerr << "EwaldgpuForce needs access to charges on GPU!" << std::endl;

  //Initialization values
  m_rcut = rcut;
  m_num_kx = num_kx;
  m_num_ky = num_ky;
  m_num_kz = num_kz;
  m_alpha = alpha;
  m_isTuned = false;

  //Compute the number of k's in k-sphere
  compute_num_k();

  set_params(m_rcut, m_num_kx, m_num_ky, m_num_kz, m_alpha);
  if(ewaldgpu_params.time_calc_steps==-1) ewaldgpu_params.time_calc_steps = determine_calc_time_steps();
}
EwaldgpuForce::~EwaldgpuForce()
{
  if(m_k)
    HANDLE_ERROR(cudaFree(m_k));
  HANDLE_ERROR(cudaFree(m_dev_k));
  if(m_forces_reci)
    HANDLE_ERROR(cudaFree(m_forces_reci));
  HANDLE_ERROR(cudaFree(m_dev_forces_reci));
  if(m_infl_factor)
    HANDLE_ERROR(cudaFree(m_infl_factor));
  HANDLE_ERROR(cudaFree(m_dev_infl_factor));
  if(m_rho_hat)
    HANDLE_ERROR(cudaFree(m_rho_hat));
  HANDLE_ERROR(cudaFree(m_dev_rho_hat));
  if(m_energy_reci)
    HANDLE_ERROR(cudaFree(m_energy_reci));
  if(m_dev_energy_reci)
    HANDLE_ERROR(cudaFree(m_dev_energy_reci));
  if(m_q_sqr)
    HANDLE_ERROR(cudaFree(m_q_sqr));
  if(m_dev_q_sqr)
    HANDLE_ERROR(cudaFree(m_dev_q_sqr));
}
void EwaldgpuForce::setup(SystemInterface &s)
{
  if (s.npart_gpu() == m_N and ewaldgpu_params.isTunedFlag) // unchanged
    {
      return;
    }

  ewaldgpu_params.isTunedFlag = ewaldgpu_params.isTuned;

  //Initialization values
  m_rcut = ewaldgpu_params.rcut;
  m_num_kx = ewaldgpu_params.num_kx;
  m_num_ky = ewaldgpu_params.num_ky;
  m_num_kz = ewaldgpu_params.num_kz;
  m_alpha = ewaldgpu_params.alpha;
  m_coulomb_prefactor = coulomb.prefactor;

  //Compute the number of k's in k-sphere
  compute_num_k();

  if (m_dev_k)
    HANDLE_ERROR(cudaFree(m_dev_k));
  HANDLE_ERROR(cudaMalloc((void**)&(m_dev_k),3*m_num_k*sizeof(real)));
  if (m_k)
    HANDLE_ERROR(cudaFreeHost(m_k));
  HANDLE_ERROR(cudaMallocHost((void**)&(m_k),3*m_num_k*sizeof(real)));
  if (m_dev_rho_hat)
    HANDLE_ERROR(cudaFree(m_dev_rho_hat));
  HANDLE_ERROR(cudaMalloc((void**)&(m_dev_rho_hat),2*m_num_k*sizeof(real)));
  if (m_rho_hat)
    HANDLE_ERROR(cudaFreeHost(m_rho_hat));
  HANDLE_ERROR(cudaMallocHost((void**)&(m_rho_hat),2*m_num_k*sizeof(real)));
  if (m_dev_infl_factor)
    HANDLE_ERROR(cudaFree(m_dev_infl_factor));
  HANDLE_ERROR(cudaMalloc((void**)&(m_dev_infl_factor),m_num_k*sizeof(real)));
  if (m_infl_factor)
    HANDLE_ERROR(cudaFreeHost(m_infl_factor));
  HANDLE_ERROR(cudaMallocHost((void**)&(m_infl_factor),m_num_k*sizeof(real)));
  if (m_dev_energy_reci)
    HANDLE_ERROR(cudaFree(m_dev_energy_reci));
  HANDLE_ERROR(cudaMalloc((void**)&(m_dev_energy_reci),sizeof(real)));
  if (m_energy_reci)
    HANDLE_ERROR(cudaFreeHost(m_energy_reci));
  HANDLE_ERROR(cudaMallocHost((void**)&(m_energy_reci),sizeof(real)));
  if (m_dev_q_sqr)
    HANDLE_ERROR(cudaFree(m_dev_q_sqr));
  HANDLE_ERROR(cudaMalloc((void**)&(m_dev_q_sqr),sizeof(real)));
  if (m_q_sqr)
    HANDLE_ERROR(cudaFreeHost(m_q_sqr));
  HANDLE_ERROR(cudaMallocHost((void**)&(m_q_sqr),sizeof(real)));

  //Resize box
  m_box_l[0] = s.box()[0];
  m_box_l[1] = s.box()[1];
  m_box_l[2] = s.box()[2];


  //Particle number
  m_N = s.npart_gpu();

  //Compute reciprocal space vectors k
  m_V=m_box_l[0]*m_box_l[1]*m_box_l[2];
  compute_k_AND_influence_factor();

  //Init GPU stream
  stream0 = (cudaStream_t *) malloc (1 * sizeof(cudaStream_t));
  start = (cudaEvent_t *) malloc (1 * sizeof(cudaEvent_t));
  stop = (cudaEvent_t *) malloc (1 * sizeof(cudaEvent_t));
  HANDLE_ERROR(cudaEventCreate(start));
  HANDLE_ERROR(cudaEventCreate(stop));
  HANDLE_ERROR(cudaStreamCreate(stream0));
  HANDLE_ERROR(cudaDeviceSynchronize());
  HANDLE_ERROR(cudaEventRecord(*start, 0));

  //Parameters
  set_params(m_rcut, m_num_kx, m_num_ky, m_num_kz, m_alpha);

  //q squared
  *m_q_sqr = 0;
  HANDLE_ERROR(cudaMemcpyAsync( m_dev_q_sqr, m_q_sqr, sizeof(real), cudaMemcpyHostToDevice));
  GPU_q_sqr(s);

  //Copy arrays on the GPU
  HANDLE_ERROR( cudaMemcpyAsync( m_dev_k, m_k, 3*m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0 ));
  HANDLE_ERROR( cudaMemcpyAsync( m_dev_infl_factor, m_infl_factor, m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
}

//Compute forces and energy
void EwaldgpuForce::computeForces(SystemInterface &s)
{
  if (coulomb.method != COULOMB_EWALD_GPU) // EWALDGPU was disabled. nobody cares about our calculations anymore
    return;

  setup(s);

  //Resize box
  m_box_l[0] = s.box()[0];
  m_box_l[1] = s.box()[1];
  m_box_l[2] = s.box()[2];

  HANDLE_ERROR( cudaMemset(m_dev_rho_hat, 0, 2*m_num_k*sizeof(real) ) );

  //Start GPU calculation
  GPU_Forces(s);
}

void EwaldgpuForce::computeEnergy(SystemInterface &s)
{
  if (coulomb.method != COULOMB_EWALD_GPU) // EWALDGPU was disabled. nobody cares about our calculations anymore
      return;

  setup(s);

  //Resize box
  m_box_l[0] = s.box()[0];
  m_box_l[1] = s.box()[1];
  m_box_l[2] = s.box()[2];
  //Set to 0
  memset(m_energy_reci,0,sizeof(real));
  //Copy arrays on the GPU
  HANDLE_ERROR( cudaMemset(m_dev_energy_reci, 0, sizeof(real)));

  //Start GPU calculation
  GPU_Energy(s);
  //Self energy
  EwaldCPU_EnergySelf();
  //Total energy
  m_energy_tot = m_energy_reci[0] +  m_energy_self;

  HANDLE_ERROR( cudaMemcpy(&(((CUDA_energy*)s.eGpu())->coulomb), &m_energy_tot,sizeof(real),cudaMemcpyHostToDevice ) );
}

//Kernel calls
void EwaldgpuForce::GPU_Forces(SystemInterface &s)
{
  //Maximum Blocks/Threads
  int maxThreadsPerBlockStructurFactor=256;
  int maxThreadsPerBlockForce=128;

  /*Kernel*/
  //Blocks, threads
  int threads;
  int blocks;
  cudaDeviceProp prop;
  HANDLE_ERROR( cudaGetDeviceProperties( &prop, 0 ) );

  /********************************************************************************************
		 	 	 	 	 	 	 	 	 	 	 	 	 	 	Structure factor / Rho_hat
  ********************************************************************************************/

  //Blocks, threads
  getNumBlocksAndThreads(m_N,  prop.sharedMemPerBlock,  maxThreadsPerBlockStructurFactor,  blocks,  threads);
  blocks=(int)ceil(sqrt(m_num_k));
  dim3 dimBlock1(threads, 1, 1);
  dim3 dimGrid1(blocks, blocks, 1);
  //Shared memory size
  int smemSize = 2 * 2 * threads * sizeof(real);

  //Loops needed in EwaldGPU_Rho_hat_MaxThreads
  int loops=m_N/(2*maxThreadsPerBlockStructurFactor);

  //Kernel call
  for(int l=0;l<loops;l++)
    {
      switch(maxThreadsPerBlockStructurFactor)
	{
	case 1024:	EwaldGPU_Rho_hat_MaxThreads<1024,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  512:	EwaldGPU_Rho_hat_MaxThreads<512,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  256:	EwaldGPU_Rho_hat_MaxThreads<256,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  128:	EwaldGPU_Rho_hat_MaxThreads<128,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  64:	EwaldGPU_Rho_hat_MaxThreads<64,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  32:	EwaldGPU_Rho_hat_MaxThreads<32,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  16:	EwaldGPU_Rho_hat_MaxThreads<16,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  8:	EwaldGPU_Rho_hat_MaxThreads<8,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  4:	EwaldGPU_Rho_hat_MaxThreads<4,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  2:	EwaldGPU_Rho_hat_MaxThreads<2,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	}
    }

  //Blocks, threads
  getNumBlocksAndThreads(m_N-2*loops*maxThreadsPerBlockStructurFactor,  prop.sharedMemPerBlock,  maxThreadsPerBlockStructurFactor,  blocks,  threads);
  blocks=(int)ceil(sqrt(m_num_k));
  dim3 dimBlock2(threads, 1, 1);
  dim3 dimGrid2(blocks, blocks, 1);

  //Kernel call
  if (isPow2(m_N-2*loops*maxThreadsPerBlockStructurFactor) && (m_N-2*loops*maxThreadsPerBlockStructurFactor != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_Rho_hat_LowThreads<1024,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_Rho_hat_LowThreads<512,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_Rho_hat_LowThreads<256,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_Rho_hat_LowThreads<128,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_Rho_hat_LowThreads<64,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_Rho_hat_LowThreads<32,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_Rho_hat_LowThreads<16,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_Rho_hat_LowThreads<8,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_Rho_hat_LowThreads<4,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_Rho_hat_LowThreads<2,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_Rho_hat_LowThreads<1,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	}
    }
  if (!isPow2(m_N-2*loops*maxThreadsPerBlockStructurFactor) && (m_N-2*loops*maxThreadsPerBlockStructurFactor != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_Rho_hat_LowThreads<1024,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_Rho_hat_LowThreads<512,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_Rho_hat_LowThreads<256,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_Rho_hat_LowThreads<128,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_Rho_hat_LowThreads<64,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_Rho_hat_LowThreads<32,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_Rho_hat_LowThreads<16,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_Rho_hat_LowThreads<8,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_Rho_hat_LowThreads<4,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_Rho_hat_LowThreads<2,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_Rho_hat_LowThreads<1,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	}
    }

  //Copy the arrays back from the GPU to the CPU
  HANDLE_ERROR(cudaMemcpy( m_rho_hat, m_dev_rho_hat,2*m_num_k*sizeof(real),cudaMemcpyDeviceToHost));

  /********************************************************************************************
																			 Forces long range
  ********************************************************************************************/

  //Blocks, threads
  getNumBlocksAndThreads(m_num_k,  prop.sharedMemPerBlock,  maxThreadsPerBlockForce,  blocks,  threads);
  blocks=(int)ceil(sqrt(m_N));
  dim3 dimBlock3(threads, 1, 1);
  dim3 dimGrid3(blocks, blocks, 1);

  //Shared memory size
  smemSize = 3 * 2 * threads * sizeof(real);

  //Loops needed in EwaldGPU_ForcesReci_MaxThreads
  loops=m_num_k/(2*maxThreadsPerBlockForce);

  //Kernel call
  for(int l=0;l<loops;l++)
    {
      switch(maxThreadsPerBlockForce)
	{
	case 1024:	EwaldGPU_ForcesReci_MaxThreads<1024,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  512:	EwaldGPU_ForcesReci_MaxThreads<512,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  256:	EwaldGPU_ForcesReci_MaxThreads<256,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  128:	EwaldGPU_ForcesReci_MaxThreads<128,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  64:	EwaldGPU_ForcesReci_MaxThreads<64,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  32:	EwaldGPU_ForcesReci_MaxThreads<32,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  16:	EwaldGPU_ForcesReci_MaxThreads<16,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  8:	EwaldGPU_ForcesReci_MaxThreads<8,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  4:	EwaldGPU_ForcesReci_MaxThreads<4,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	case  2:	EwaldGPU_ForcesReci_MaxThreads<2,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,l);cuda_check_error(dimBlock3, dimGrid3, 0, __FILE__, __LINE__);break;
	}
    }

  //Blocks, threads
  getNumBlocksAndThreads(m_num_k-2*loops*maxThreadsPerBlockForce,  prop.sharedMemPerBlock,  maxThreadsPerBlockForce,  blocks,  threads);
  blocks=(int)ceil(sqrt(m_N));
  dim3 dimBlock4(threads, 1, 1);
  dim3 dimGrid4(blocks, blocks, 1);

  //Shared memory size
  smemSize = 3 * 2 * threads * sizeof(real);

  //Kernel call
  if (isPow2(m_num_k-2*loops*maxThreadsPerBlockForce) && (m_num_k-2*loops*maxThreadsPerBlockForce != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_ForcesReci_LowThreads<1024,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_ForcesReci_LowThreads<512,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_ForcesReci_LowThreads<256,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_ForcesReci_LowThreads<128,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_ForcesReci_LowThreads<64,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_ForcesReci_LowThreads<32,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_ForcesReci_LowThreads<16,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_ForcesReci_LowThreads<8,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_ForcesReci_LowThreads<4,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_ForcesReci_LowThreads<2,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_ForcesReci_LowThreads<1,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	}
    }
  if (!isPow2(m_num_k-2*loops*maxThreadsPerBlockForce) && (m_num_k-2*loops*maxThreadsPerBlockForce != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_ForcesReci_LowThreads<1024,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_ForcesReci_LowThreads<512,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_ForcesReci_LowThreads<256,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_ForcesReci_LowThreads<128,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_ForcesReci_LowThreads<64,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_ForcesReci_LowThreads<32,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_ForcesReci_LowThreads<16,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_ForcesReci_LowThreads<8,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_ForcesReci_LowThreads<4,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_ForcesReci_LowThreads<2,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_ForcesReci_LowThreads<1,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,s.rGpuBegin(),s.qGpuBegin(),m_dev_rho_hat,m_dev_infl_factor,s.fGpuBegin(),m_N,m_num_k,m_V,m_coulomb_prefactor,maxThreadsPerBlockForce,loops); cuda_check_error(dimBlock4, dimGrid4, 0, __FILE__, __LINE__);break;
	}
    }
}
void EwaldgpuForce::GPU_Energy(SystemInterface &s)
{
  //Maximum Blocks/Threads
  int maxThreadsPerBlockEnergie=128;

  /*Kernel*/
  //Blocks, threads
  int threads;
  int blocks;
  cudaDeviceProp prop;
  HANDLE_ERROR( cudaGetDeviceProperties( &prop, 0 ) );

  /********************************************************************************************
																					 Energy
  ********************************************************************************************/

  //Blocks, threads
  getNumBlocksAndThreads(m_num_k,  prop.sharedMemPerBlock,  maxThreadsPerBlockEnergie,  blocks,  threads);
  blocks=1;
  dim3 dimBlock1(threads, 1, 1);
  dim3 dimGrid1(blocks, 1, 1);

  //Shared memory size
  int smemSize = 1 * 2 * threads * sizeof(real);

  //Loops needed in EwaldGPU_ForcesReci_MaxThreads
  int loops=m_num_k/(2*maxThreadsPerBlockEnergie);

  //Kernel call
  for(int l=0;l<loops;l++)
    {
      switch(maxThreadsPerBlockEnergie)
	{
	case 1024:	EwaldGPU_EnergyReci_MaxThreads<1024,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  512:	EwaldGPU_EnergyReci_MaxThreads<512,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  256:	EwaldGPU_EnergyReci_MaxThreads<256,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  128:	EwaldGPU_EnergyReci_MaxThreads<128,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  64:	EwaldGPU_EnergyReci_MaxThreads<64,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  32:	EwaldGPU_EnergyReci_MaxThreads<32,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  16:	EwaldGPU_EnergyReci_MaxThreads<16,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  8:	EwaldGPU_EnergyReci_MaxThreads<8,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  4:	EwaldGPU_EnergyReci_MaxThreads<4,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  2:	EwaldGPU_EnergyReci_MaxThreads<2,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l,m_coulomb_prefactor);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	}
    }

  //Blocks, threads
  getNumBlocksAndThreads(m_num_k-2*loops*maxThreadsPerBlockEnergie,  prop.sharedMemPerBlock,  maxThreadsPerBlockEnergie,  blocks,  threads);
  blocks=1;
  dim3 dimBlock2(threads, 1, 1);
  dim3 dimGrid2(blocks, 1, 1);

  //Shared memory size
  smemSize = 1 * 2 * threads * sizeof(real);

  //Kernel call
  if (isPow2(m_num_k-2*loops*maxThreadsPerBlockEnergie) && (m_num_k-2*loops*maxThreadsPerBlockEnergie != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_EnergyReci_LowThreads<1024,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_EnergyReci_LowThreads<512,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_EnergyReci_LowThreads<256,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_EnergyReci_LowThreads<128,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_EnergyReci_LowThreads<64,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_EnergyReci_LowThreads<32,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_EnergyReci_LowThreads<16,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_EnergyReci_LowThreads<8,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_EnergyReci_LowThreads<4,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_EnergyReci_LowThreads<2,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_EnergyReci_LowThreads<1,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	}
    }
  if (!isPow2(m_num_k-2*loops*maxThreadsPerBlockEnergie) && (m_num_k-2*loops*maxThreadsPerBlockEnergie != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_EnergyReci_LowThreads<1024,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_EnergyReci_LowThreads<512,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_EnergyReci_LowThreads<256,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_EnergyReci_LowThreads<128,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_EnergyReci_LowThreads<64,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_EnergyReci_LowThreads<32,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_EnergyReci_LowThreads<16,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_EnergyReci_LowThreads<8,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_EnergyReci_LowThreads<4,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_EnergyReci_LowThreads<2,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_EnergyReci_LowThreads<1,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,&(((CUDA_energy*)s.eGpu())->coulomb),m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops,m_coulomb_prefactor); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	}
    }

  //Copy the values back from the GPU to the CPU
  HANDLE_ERROR( cudaMemcpy( m_energy_reci, &(((CUDA_energy*)s.eGpu())->coulomb),sizeof(real),cudaMemcpyDeviceToHost ) );
}

struct square {
  __host__ __device__ 
  real operator()(const real &x) const {
    return x*x;
  }
};

void EwaldgpuForce::GPU_q_sqr(SystemInterface &s)
{
  //Maximum Blocks/Threads
  int maxThreadsPerBlock_q_sqr=128;

  /*Kernel*/
  //Blocks, threads
  int threads;
  int blocks;
  cudaDeviceProp prop;
  HANDLE_ERROR( cudaGetDeviceProperties( &prop, 0 ) );

  /********************************************************************************************
																					 q squared
  ********************************************************************************************/

  //Blocks, threads
  getNumBlocksAndThreads(m_N,  prop.sharedMemPerBlock,  maxThreadsPerBlock_q_sqr,  blocks,  threads);
  blocks=1;
  dim3 dimBlock1(threads, 1, 1);
  dim3 dimGrid1(blocks, 1, 1);

  //Shared memory size
  int smemSize = 1 * 2 * threads * sizeof(real);

  //Loops needed in EwaldGPU_ForcesReci_MaxThreads
  int loops=m_N/(2*maxThreadsPerBlock_q_sqr);

  //Kernel call
  for(int l=0;l<loops;l++)
    {
      switch(maxThreadsPerBlock_q_sqr)
	{
	case 1024:	EwaldGPU_q_sqr_MaxThreads<1024,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  512:	EwaldGPU_q_sqr_MaxThreads<512,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  256:	EwaldGPU_q_sqr_MaxThreads<256,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  128:	EwaldGPU_q_sqr_MaxThreads<128,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  64:	EwaldGPU_q_sqr_MaxThreads<64,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  32:	EwaldGPU_q_sqr_MaxThreads<32,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  16:	EwaldGPU_q_sqr_MaxThreads<16,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  8:	EwaldGPU_q_sqr_MaxThreads<8,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  4:	EwaldGPU_q_sqr_MaxThreads<4,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	case  2:	EwaldGPU_q_sqr_MaxThreads<2,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,l);cuda_check_error(dimBlock1, dimGrid1, 0, __FILE__, __LINE__);break;
	}
    }

  //Blocks, threads
  getNumBlocksAndThreads(m_N-2*loops*maxThreadsPerBlock_q_sqr,  prop.sharedMemPerBlock,  maxThreadsPerBlock_q_sqr,  blocks,  threads);
  blocks=1;
  dim3 dimBlock2(threads, 1, 1);
  dim3 dimGrid2(blocks, 1, 1);

  //Shared memory size
  smemSize = 1 * 2 * threads * sizeof(real);

  //Kernel call
  if (isPow2(m_N-2*loops*maxThreadsPerBlock_q_sqr) && (m_N-2*loops*maxThreadsPerBlock_q_sqr != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_q_sqr_LowThreads<1024,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_q_sqr_LowThreads<512,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_q_sqr_LowThreads<256,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_q_sqr_LowThreads<128,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_q_sqr_LowThreads<64,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_q_sqr_LowThreads<32,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_q_sqr_LowThreads<16,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_q_sqr_LowThreads<8,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_q_sqr_LowThreads<4,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_q_sqr_LowThreads<2,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_q_sqr_LowThreads<1,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	}
    }
  if (!isPow2(m_N-2*loops*maxThreadsPerBlock_q_sqr) && (m_N-2*loops*maxThreadsPerBlock_q_sqr != 0))
    {
      switch (threads)
	{
	case 1024: EwaldGPU_q_sqr_LowThreads<1024,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  512: EwaldGPU_q_sqr_LowThreads<512,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  256: EwaldGPU_q_sqr_LowThreads<256,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case  128: EwaldGPU_q_sqr_LowThreads<128,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   64: EwaldGPU_q_sqr_LowThreads<64,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   32: EwaldGPU_q_sqr_LowThreads<32,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case   16: EwaldGPU_q_sqr_LowThreads<16,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    8: EwaldGPU_q_sqr_LowThreads<8,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    4: EwaldGPU_q_sqr_LowThreads<4,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    2: EwaldGPU_q_sqr_LowThreads<2,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	case    1: EwaldGPU_q_sqr_LowThreads<1,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(s.qGpuBegin(),m_dev_q_sqr,m_N,maxThreadsPerBlock_q_sqr,loops); cuda_check_error(dimBlock2, dimGrid2, 0, __FILE__, __LINE__);break;
	}
    }

  //Copy the values back from the GPU to the CPU
  HANDLE_ERROR( cudaMemcpy( m_q_sqr, m_dev_q_sqr,sizeof(real),cudaMemcpyDeviceToHost) );
}
void cuda_check_error(const dim3 &block, const dim3 &grid, const char *function, const char *file, unsigned int line)
{
  _err=cudaGetLastError();
  if (_err!=cudaSuccess)
    {
      fprintf(stderr, "%d: error \"%s\" calling %s with dim %d %d %d, grid %d %d %d in %s:%u\n", this_node, cudaGetErrorString(_err), function, block.x, block.y, block.z, grid.x, grid.y, grid.z,file, line);
      exit(EXIT_FAILURE);
    }
}

//Compute blocks and threads
int  EwaldgpuForce::nextPow2(int x)
{
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}
bool EwaldgpuForce::isPow2(int x)
{
  if(x==1) return false;
  else return ((x&(x-1))==0);
}
void EwaldgpuForce::getNumBlocksAndThreads(int Size, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
  //Get device capability, to avoid block/grid size exceed the upbound
  cudaDeviceProp prop;
  int device;

  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&prop,device));

  threads = (Size < maxThreads*2) ? nextPow2((Size + 1)/ 2) : maxThreads;
}

#endif
