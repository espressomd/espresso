#include "config.hpp"
#ifdef EWALD_GPU

#include <mpi.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include "interaction_data.hpp"
#include "grid.hpp"
#include "tuning.hpp"
#include "integrate.hpp"
#include "EwaldgpuForce.hpp"
typedef ewaldgpu_real real;
Ewaldgpu_params ewaldgpu_params;

//Error handler
static void HandleError(cudaError_t err,const char *file,int line)
{
    if (err != cudaSuccess)
    {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define HANDLE_NULL( a ) {if (a == NULL) {printf( "Host memory failed in %s at line %d\n", __FILE__, __LINE__ ); exit( EXIT_FAILURE );}}

/*Here MaxThreads/LowThreads means that if for example 2000 Particles/k-Vectors are given and the GPU uses 256 Threads the first 3*(2*256)=1536 will be
reduced by the MaxThread function and the remaining 464 Particles will be reduced by the LowThread function*/
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
	{
		//Init
		i = tid;
		sdata[tid] = 0;
		sdata[tid+2*blockSize] = 0;

		//REDUCTION
		while (i < mTPB)
		{
			kr = k[blockId]*r_i[i+2*l*mTPB]+k[blockId+num_k]*r_i[i+2*l*mTPB+N]+k[blockId+2*num_k]*r_i[i+2*l*mTPB+2*N];
			sincos(kr,sin_ptr,cos_ptr);
			factor = q_i[i+2*l*mTPB];

			sdata[tid]      						+=  factor*cos_kr;
			sdata[tid+2*blockSize]      += -factor*sin_kr;
			//BECAUSE nIsPow2=True
			kr = k[blockId]*r_i[i+2*l*mTPB+blockSize]+k[blockId+num_k]*r_i[i+2*l*mTPB+blockSize+N]+k[blockId+2*num_k]*r_i[i+2*l*mTPB+blockSize+2*N];
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

	//REDUCTION
	while (i < N-2*l*mTPB)
	{
		kr = k[blockId]*r_i[i+2*l*mTPB]+k[blockId+num_k]*r_i[i+2*l*mTPB+N]+k[blockId+2*num_k]*r_i[i+2*l*mTPB+2*N];
		sincos(kr,sin_ptr,cos_ptr);
		factor = q_i[i+2*l*mTPB];

		sdata[tid]      						+=  factor*cos_kr;
		sdata[tid+2*blockSize]      += -factor*sin_kr;
		if (nIsPow2 || i + blockSize < N-2*l*mTPB)
		{
		kr = k[blockId]*r_i[i+2*l*mTPB+blockSize]+k[blockId+num_k]*r_i[i+2*l*mTPB+blockSize+N]+k[blockId+2*num_k]*r_i[i+2*l*mTPB+blockSize+2*N];
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
__global__ void EwaldGPU_ForcesReci_MaxThreads(real *k,real *r_i, real *q_i, real *rho_hat, real *infl_factor, real *forces_reci, int N, int num_k, real V,int maxThreadsPerBlock,int loops)
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
	{
		//Init
		i = tid;
		sdata[tid] 						 = 0;
		sdata[tid+2*blockSize] = 0;
		sdata[tid+4*blockSize] = 0;

		//REDUCTION
		while (i < mTPB)
		{
			kr = k[i+2*l*mTPB]*r_i[blockId] + k[i+2*l*mTPB+num_k]*r_i[blockId+N] + k[i+2*l*mTPB+2*num_k]*r_i[blockId+2*N];
			sincos(kr,sin_ptr,cos_ptr);
			factor = infl_factor[i+2*l*mTPB] * q_i[blockId];

			sdata[tid]      			 += factor * (k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB]         * sin_kr + k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB+num_k]         * cos_kr);
			sdata[tid+2*blockSize] += factor * (k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB]   * sin_kr + k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB+num_k]   * cos_kr);
			sdata[tid+4*blockSize] += factor * (k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB] * sin_kr + k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB+num_k] * cos_kr);
			//BECAUSE nIsPow2=True
			kr = k[i+2*l*mTPB+blockSize]*r_i[blockId] + k[i+2*l*mTPB+num_k+blockSize]*r_i[blockId+N] + k[i+2*l*mTPB+2*num_k+blockSize]*r_i[blockId+2*N];
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
			forces_reci[blockId]     += sdata[0];
			forces_reci[blockId+N]   += sdata[2*blockSize];
			forces_reci[blockId+2*N] += sdata[4*blockSize];
		}
		__syncthreads();
	}
}
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_ForcesReci_LowThreads(real *k,real *r_i, real *q_i, real *rho_hat, real *infl_factor, real *forces_reci, int N, int num_k, real V,int maxThreadsPerBlock,int elapsedLoops)
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
		kr = k[i+2*l*mTPB]*r_i[blockId] + k[i+2*l*mTPB+num_k]*r_i[blockId+N] + k[i+2*l*mTPB+2*num_k]*r_i[blockId+2*N];
		sincos(kr,sin_ptr,cos_ptr);
		factor = infl_factor[i+2*l*mTPB] * q_i[blockId];

		//REDUCTION
		sdata[tid]      			 += factor * (k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB]         * sin_kr + k[i+2*l*mTPB]*rho_hat[i+2*l*mTPB+num_k]         * cos_kr);
		sdata[tid+2*blockSize] += factor * (k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB]   * sin_kr + k[i+2*l*mTPB+num_k]*rho_hat[i+2*l*mTPB+num_k]   * cos_kr);
		sdata[tid+4*blockSize] += factor * (k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB] * sin_kr + k[i+2*l*mTPB+2*num_k]*rho_hat[i+2*l*mTPB+num_k] * cos_kr);
		if (nIsPow2 || i + blockSize < num_k-2*l*mTPB)
		{
			kr = k[i+2*l*mTPB+blockSize]*r_i[blockId] + k[i+2*l*mTPB+num_k+blockSize]*r_i[blockId+N] + k[i+2*l*mTPB+2*num_k+blockSize]*r_i[blockId+2*N];
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
		forces_reci[blockId]     += sdata[0];
		forces_reci[blockId+N]   += sdata[2*blockSize];
		forces_reci[blockId+2*N] += sdata[4*blockSize];
	}
}
//Energy in reciprocal space
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_EnergyReci_MaxThreads(real *rho_hat, real *infl_factor,real *energy_reci, int N, int num_k, real V,int maxThreadsPerBlock,int loops)
{
	//Variables
	extern __shared__ real sdata[];
	int tid = threadIdx.x;
	int i = tid;
	int gridSize = blockSize*2*gridDim.x;
	int mTPB=maxThreadsPerBlock;
	real factor; // Factor to multiply with

	int l = loops;//TODO
	{
		//iNIT
		i = tid;
		sdata[tid] 						 = 0;

		//REDUCTION
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
			energy_reci[0] += sdata[0];
		}
		__syncthreads();
	}
}
template <int blockSize, bool nIsPow2>
__global__ void EwaldGPU_EnergyReci_LowThreads(real *rho_hat, real *infl_factor,real *energy_reci, int N, int num_k, real V,int maxThreadsPerBlock,int elapsedLoops)
{
	//Variables
	extern __shared__ real sdata[];
	int tid = threadIdx.x;
	int i = tid;
	int gridSize = blockSize*2*gridDim.x;
	int mTPB=maxThreadsPerBlock;
	int l=elapsedLoops;
	real factor;

	//iNIT
	i = tid;
	sdata[tid] 						 = 0;

	//REDUCTION
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
		energy_reci[0] += sdata[0];
	}
}

//EWALDGPUFORCE
EwaldgpuForce::EwaldgpuForce(double rcut, int num_kx, int num_ky, int num_kz, double alpha)
{
	//Initialization values
	m_rcut = rcut;
	m_num_kx = num_kx;
	m_num_ky = num_ky;
	m_num_kz = num_kz;
	m_alpha = alpha;

	m_coulomb_prefactor = coulomb.prefactor;
	initialized = false;

	//Compute the number of k's in k-sphere
	compute_num_k();

	//COULOMB METHOD
	coulomb.method = COULOMB_EWALD_GPU;
	ewaldgpu_set_params(m_rcut, m_num_kx, m_num_ky, m_num_kz, m_alpha);
	if(ewaldgpu_params.time_calc_steps==0) ewaldgpu_params.time_calc_steps = ewaldgpu_count_charged_particles();
	mpi_bcast_coulomb_params();
}
EwaldgpuForce::~EwaldgpuForce()
{
	//Free memory
	if (initialized) destroy();

	//RESET GPU
	HANDLE_ERROR(cudaDeviceReset());
	HANDLE_ERROR(cudaFree(0));
}
void EwaldgpuForce::init(SystemInterface &systemInterface)
{
	//System interface
	System = &systemInterface;

	// If number of particles has changed reinitialize
	if (m_N == systemInterface.npart() and isTuned and ewaldgpu_params.isTuned) return;
	isTuned = ewaldgpu_params.isTuned;

	//Initialization values
	m_rcut = ewaldgpu_params.rcut;
	m_num_kx = ewaldgpu_params.num_kx;
	m_num_ky = ewaldgpu_params.num_ky;
	m_num_kz = ewaldgpu_params.num_kz;
	m_alpha = ewaldgpu_params.alpha;
	//Compute the number of k's in k-sphere
	compute_num_k();

	//MEMORY CPU
	//Free memory and destroy events and streams if already initialized
  if(initialized) destroy();
	initialized=true;

	//Reserve memory for force vector
	m_N = systemInterface.npart();
	F.reserve(systemInterface.npart());

	//Allocate memory on the CPU
	HANDLE_ERROR(cudaMallocHost((void**)&(m_k),3*m_num_k*sizeof(real)));
	HANDLE_ERROR(cudaMallocHost((void**)&(m_rho_hat),2*m_num_k*sizeof(real)));
	HANDLE_ERROR(cudaMallocHost((void**)&(m_q_i),m_N*sizeof(real)));
	HANDLE_ERROR(cudaMallocHost((void**)&(m_r_i),3*m_N*sizeof(real)));
	HANDLE_ERROR(cudaMallocHost((void**)&(m_infl_factor),m_num_k*sizeof(real)));
	HANDLE_ERROR(cudaMallocHost((void**)&(m_forces_reci),3*m_N*sizeof(real)));
	HANDLE_ERROR(cudaMallocHost((void**)&(m_energy_reci),sizeof(real)));
	m_energy_self = (real*)malloc(sizeof(real));


	//Transform the positions and currents from Vector3d-format to array format needed on GPU
	Transform_Vector3d_to_Array(systemInterface);

	//COMPUTE REZIPROCAL K's
	m_V=m_box_l[0]*m_box_l[1]*m_box_l[2];
	compute_q_sqare();
	compute_k_AND_influence_factor();


	//INIT GPU STREAM
	stream0 = (cudaStream_t *) malloc (1 * sizeof(cudaStream_t));
	start = (cudaEvent_t *) malloc (1 * sizeof(cudaEvent_t));
	stop = (cudaEvent_t *) malloc (1 * sizeof(cudaEvent_t));
	HANDLE_ERROR(cudaEventCreate(&(*start)));
	HANDLE_ERROR(cudaEventCreate(&(*stop)));
	HANDLE_ERROR(cudaStreamCreate(&(*stream0)));
	HANDLE_ERROR(cudaDeviceSynchronize());
	cudaEventRecord(*start, 0);

  //MEMORY GPU
  //Allocate memory on the GPU
	HANDLE_ERROR(cudaMalloc((void**)&(m_dev_k),3*m_num_k*sizeof(real)));
	HANDLE_ERROR(cudaMalloc((void**)&(m_dev_rho_hat),2*m_num_k*sizeof(real)));
	HANDLE_ERROR(cudaMalloc((void**)&(m_dev_q_i),m_N*sizeof(real)));
	HANDLE_ERROR(cudaMalloc((void**)&(m_dev_r_i),3*m_N*sizeof(real)));
	HANDLE_ERROR(cudaMalloc((void**)&(m_dev_infl_factor),m_num_k*sizeof(real)));
	HANDLE_ERROR(cudaMalloc((void**)&(m_dev_forces_reci),3*m_N*sizeof(real)));
	HANDLE_ERROR(cudaMalloc((void**)&(m_dev_energy_reci),sizeof(real)));

	ewaldgpu_set_params(m_rcut, m_num_kx, m_num_ky, m_num_kz, m_alpha);
}
void EwaldgpuForce::run(SystemInterface &systemInterface)
{
	if (coulomb.method != COULOMB_EWALD_GPU)
	{
		printf("Error: It is currently not supported to disable forces using the EspressoSystemInterface.\n");
		exit(EXIT_FAILURE);
	}
	//Transform the positions and currents from Vector3d-format to array format needed on GPU
	Transform_Vector3d_to_Array(systemInterface);
	//Set to NULL
	memset(m_rho_hat,0,2*m_num_k*sizeof(real));
	memset(m_forces_reci,0,3*m_N*sizeof(real));
	memset(m_energy_reci,0,sizeof(real));
	//Copy arrays on the GPU
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_k, m_k, 3*m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0 ));
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_q_i, m_q_i, m_N*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_r_i, m_r_i, 3*m_N*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_rho_hat, m_rho_hat, 2*m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_infl_factor, m_infl_factor, m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_forces_reci, m_forces_reci, 3*m_N*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_energy_reci, m_energy_reci, sizeof(real), cudaMemcpyHostToDevice, 0  ) );

	//Start GPU calculation
	runGPU_Forces();
}
bool EwaldgpuForce::isReady()
{
  //Wait for GPU stream end
 	cudaEventRecord(*stop, 0);
	memset(m_energy_self,0,sizeof(real));
	EwaldCPU_EnergySelf();
  HANDLE_ERROR(cudaStreamSynchronize(*stream0));
  HANDLE_ERROR(cudaDeviceSynchronize());//TODO

	//Transforms force result to Eigen::Vector3d used in TCL
	Transform_ForceArray_to_Vector3d();

	//Output();

	return true;
}
void EwaldgpuForce::runEnergies(SystemInterface &systemInterface)
{
	if (coulomb.method != COULOMB_EWALD_GPU)
	{
		printf("Error: It is currently not supported to disable forces using the EspressoSystemInterface.\n");
		exit(EXIT_FAILURE);
	}
	if (m_N != systemInterface.npart())
	{
		printf("Error: number of particles changed between init (%d) and run (%d).\n", m_N, systemInterface.npart());
		exit(EXIT_FAILURE);
	}

	//Transform the positions and currents from Vector3d-format to array format needed on GPU
	Transform_Vector3d_to_Array(systemInterface);
	//Set to NULL

	memset(m_energy_reci,0,sizeof(real));
	//Copy arrays on the GPU
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_k, m_k, 3*m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0 ));
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_q_i, m_q_i, m_N*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_r_i, m_r_i, 3*m_N*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_rho_hat, m_rho_hat, 2*m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_infl_factor, m_infl_factor, m_num_k*sizeof(real), cudaMemcpyHostToDevice, 0  ) );
	HANDLE_ERROR( cudaMemcpyAsync( m_dev_energy_reci, m_energy_reci, sizeof(real), cudaMemcpyHostToDevice, 0  ) );

	//Start GPU calculation
	runGPU_Energy();
}
bool EwaldgpuForce::isReadyEnergies()
{
	  //Wait for GPU stream end
	 	cudaEventRecord(*stop, 0);
		memset(m_energy_self,0,sizeof(real));
		EwaldCPU_EnergySelf();
	  HANDLE_ERROR(cudaStreamSynchronize(*stream0));
	  HANDLE_ERROR(cudaDeviceSynchronize());//TODO

		//Transforms force result to Eigen::Vector3d used in TCL
		Transform_ForceArray_to_Vector3d();

		//Total energy
		energy.coulomb[0] = m_energy_reci[0] +  m_energy_self[0];

		//Output();

		return true;
}
void EwaldgpuForce::destroy()
{
	//FREE MEMORY
	//Free memory on the GPU side
	HANDLE_ERROR(cudaFree(m_dev_q_i));
	HANDLE_ERROR(cudaFree(m_dev_r_i));
	HANDLE_ERROR(cudaFree(m_dev_forces_reci));
	HANDLE_ERROR(cudaFree(m_dev_energy_reci));
	HANDLE_ERROR(cudaFree(m_dev_k));
	HANDLE_ERROR(cudaFree(m_dev_rho_hat));
	HANDLE_ERROR(cudaFree(m_dev_infl_factor));
  //Free memory streams
	HANDLE_ERROR(cudaEventDestroy(*start));
	HANDLE_ERROR(cudaEventDestroy(*stop));
	HANDLE_ERROR(cudaStreamDestroy(*stream0));
	//Free memory on the CPU side
	HANDLE_ERROR(cudaFreeHost(m_q_i));
	HANDLE_ERROR(cudaFreeHost(m_r_i));
	HANDLE_ERROR(cudaFreeHost(m_forces_reci));
	HANDLE_ERROR(cudaFreeHost(m_energy_reci));
	HANDLE_ERROR(cudaFreeHost(m_k));
	HANDLE_ERROR(cudaFreeHost(m_rho_hat));
	HANDLE_ERROR(cudaFreeHost(m_infl_factor));

	initialized=false;
}

//RUN GPU PART
void EwaldgpuForce::runGPU_Forces()
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
			case 1024:	EwaldGPU_Rho_hat_MaxThreads<1024,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  512:	EwaldGPU_Rho_hat_MaxThreads<512,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  256:	EwaldGPU_Rho_hat_MaxThreads<256,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  128:	EwaldGPU_Rho_hat_MaxThreads<128,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  64:	EwaldGPU_Rho_hat_MaxThreads<64,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  32:	EwaldGPU_Rho_hat_MaxThreads<32,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  16:	EwaldGPU_Rho_hat_MaxThreads<16,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  8:	EwaldGPU_Rho_hat_MaxThreads<8,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  4:	EwaldGPU_Rho_hat_MaxThreads<4,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
			case  2:	EwaldGPU_Rho_hat_MaxThreads<2,true><<<dimGrid1, dimBlock1, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,l);break;
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
			case 1024: EwaldGPU_Rho_hat_LowThreads<1024,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case  512: EwaldGPU_Rho_hat_LowThreads<512,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case  256: EwaldGPU_Rho_hat_LowThreads<256,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case  128: EwaldGPU_Rho_hat_LowThreads<128,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case   64: EwaldGPU_Rho_hat_LowThreads<64,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case   32: EwaldGPU_Rho_hat_LowThreads<32,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case   16: EwaldGPU_Rho_hat_LowThreads<16,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    8: EwaldGPU_Rho_hat_LowThreads<8,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    4: EwaldGPU_Rho_hat_LowThreads<4,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    2: EwaldGPU_Rho_hat_LowThreads<2,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    1: EwaldGPU_Rho_hat_LowThreads<1,true><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
		}
	}
	if (!isPow2(m_N-2*loops*maxThreadsPerBlockStructurFactor) && (m_N-2*loops*maxThreadsPerBlockStructurFactor != 0))
	{
		switch (threads)
		{
			case 1024: EwaldGPU_Rho_hat_LowThreads<1024,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case  512: EwaldGPU_Rho_hat_LowThreads<512,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case  256: EwaldGPU_Rho_hat_LowThreads<256,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case  128: EwaldGPU_Rho_hat_LowThreads<128,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case   64: EwaldGPU_Rho_hat_LowThreads<64,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case   32: EwaldGPU_Rho_hat_LowThreads<32,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case   16: EwaldGPU_Rho_hat_LowThreads<16,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    8: EwaldGPU_Rho_hat_LowThreads<8,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    4: EwaldGPU_Rho_hat_LowThreads<4,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    2: EwaldGPU_Rho_hat_LowThreads<2,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
			case    1: EwaldGPU_Rho_hat_LowThreads<1,false><<<dimGrid2, dimBlock2, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_N,m_num_k,maxThreadsPerBlockStructurFactor,loops); break;
		}
	}
	//Copy the arrays back from the GPU to the CPU
	HANDLE_ERROR( cudaMemcpyAsync( m_rho_hat, m_dev_rho_hat,2*m_num_k*sizeof(real),cudaMemcpyDeviceToHost, *stream0 ) );

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
			case 1024:	EwaldGPU_ForcesReci_MaxThreads<1024,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  512:	EwaldGPU_ForcesReci_MaxThreads<512,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  256:	EwaldGPU_ForcesReci_MaxThreads<256,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  128:	EwaldGPU_ForcesReci_MaxThreads<128,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  64:	EwaldGPU_ForcesReci_MaxThreads<64,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  32:	EwaldGPU_ForcesReci_MaxThreads<32,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  16:	EwaldGPU_ForcesReci_MaxThreads<16,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  8:	EwaldGPU_ForcesReci_MaxThreads<8,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  4:	EwaldGPU_ForcesReci_MaxThreads<4,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
			case  2:	EwaldGPU_ForcesReci_MaxThreads<2,true><<<dimGrid3, dimBlock3, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,l);break;
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
			case 1024: EwaldGPU_ForcesReci_LowThreads<1024,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case  512: EwaldGPU_ForcesReci_LowThreads<512,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case  256: EwaldGPU_ForcesReci_LowThreads<256,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case  128: EwaldGPU_ForcesReci_LowThreads<128,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case   64: EwaldGPU_ForcesReci_LowThreads<64,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case   32: EwaldGPU_ForcesReci_LowThreads<32,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case   16: EwaldGPU_ForcesReci_LowThreads<16,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    8: EwaldGPU_ForcesReci_LowThreads<8,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    4: EwaldGPU_ForcesReci_LowThreads<4,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    2: EwaldGPU_ForcesReci_LowThreads<2,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    1: EwaldGPU_ForcesReci_LowThreads<1,true><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
		}
	}
	if (!isPow2(m_num_k-2*loops*maxThreadsPerBlockForce) && (m_num_k-2*loops*maxThreadsPerBlockForce != 0))
	{
		switch (threads)
		{
			case 1024: EwaldGPU_ForcesReci_LowThreads<1024,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case  512: EwaldGPU_ForcesReci_LowThreads<512,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case  256: EwaldGPU_ForcesReci_LowThreads<256,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case  128: EwaldGPU_ForcesReci_LowThreads<128,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case   64: EwaldGPU_ForcesReci_LowThreads<64,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case   32: EwaldGPU_ForcesReci_LowThreads<32,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case   16: EwaldGPU_ForcesReci_LowThreads<16,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    8: EwaldGPU_ForcesReci_LowThreads<8,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    4: EwaldGPU_ForcesReci_LowThreads<4,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    2: EwaldGPU_ForcesReci_LowThreads<2,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
			case    1: EwaldGPU_ForcesReci_LowThreads<1,false><<<dimGrid4, dimBlock4, smemSize, *stream0>>>(m_dev_k,m_dev_r_i,m_dev_q_i,m_dev_rho_hat,m_dev_infl_factor,m_dev_forces_reci,m_N,m_num_k,m_V,maxThreadsPerBlockForce,loops); break;
		}
	}
	//Copy the arrays back from the GPU to the CPU
	HANDLE_ERROR( cudaMemcpyAsync( m_forces_reci, m_dev_forces_reci,3*m_N*sizeof(real),cudaMemcpyDeviceToHost, *stream0 ) );
}
void EwaldgpuForce::runGPU_Energy()
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
	dim3 dimBlock5(threads, 1, 1);
	dim3 dimGrid5(blocks, 1, 1);

	//Shared memory size
	int smemSize = 1 * 2 * threads * sizeof(real);

	//Loops needed in EwaldGPU_ForcesReci_MaxThreads
	int loops=m_num_k/(2*maxThreadsPerBlockEnergie);

	//Kernel call
	for(int l=0;l<loops;l++)
	{
		switch(maxThreadsPerBlockEnergie)
		{
			case 1024:	EwaldGPU_EnergyReci_MaxThreads<1024,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  512:	EwaldGPU_EnergyReci_MaxThreads<512,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  256:	EwaldGPU_EnergyReci_MaxThreads<256,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  128:	EwaldGPU_EnergyReci_MaxThreads<128,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  64:	EwaldGPU_EnergyReci_MaxThreads<64,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  32:	EwaldGPU_EnergyReci_MaxThreads<32,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  16:	EwaldGPU_EnergyReci_MaxThreads<16,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  8:	EwaldGPU_EnergyReci_MaxThreads<8,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  4:	EwaldGPU_EnergyReci_MaxThreads<4,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
			case  2:	EwaldGPU_EnergyReci_MaxThreads<2,true><<<dimGrid5, dimBlock5, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,l);break;
		}
	}

	//Blocks, threads
	getNumBlocksAndThreads(m_num_k-2*loops*maxThreadsPerBlockEnergie,  prop.sharedMemPerBlock,  maxThreadsPerBlockEnergie,  blocks,  threads);
	blocks=1;
	dim3 dimBlock6(threads, 1, 1);
	dim3 dimGrid6(blocks, 1, 1);

	//Shared memory size
	smemSize = 1 * 2 * threads * sizeof(real);

	//Kernel call
	if (isPow2(m_num_k-2*loops*maxThreadsPerBlockEnergie) && (m_num_k-2*loops*maxThreadsPerBlockEnergie != 0))
	{
		switch (threads)
		{
			case 1024: EwaldGPU_EnergyReci_LowThreads<1024,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case  512: EwaldGPU_EnergyReci_LowThreads<512,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case  256: EwaldGPU_EnergyReci_LowThreads<256,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case  128: EwaldGPU_EnergyReci_LowThreads<128,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case   64: EwaldGPU_EnergyReci_LowThreads<64,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case   32: EwaldGPU_EnergyReci_LowThreads<32,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case   16: EwaldGPU_EnergyReci_LowThreads<16,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    8: EwaldGPU_EnergyReci_LowThreads<8,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    4: EwaldGPU_EnergyReci_LowThreads<4,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    2: EwaldGPU_EnergyReci_LowThreads<2,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    1: EwaldGPU_EnergyReci_LowThreads<1,true><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
		}
	}
	if (!isPow2(m_num_k-2*loops*maxThreadsPerBlockEnergie) && (m_num_k-2*loops*maxThreadsPerBlockEnergie != 0))
	{
		switch (threads)
		{
			case 1024: EwaldGPU_EnergyReci_LowThreads<1024,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case  512: EwaldGPU_EnergyReci_LowThreads<512,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case  256: EwaldGPU_EnergyReci_LowThreads<256,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case  128: EwaldGPU_EnergyReci_LowThreads<128,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case   64: EwaldGPU_EnergyReci_LowThreads<64,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case   32: EwaldGPU_EnergyReci_LowThreads<32,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case   16: EwaldGPU_EnergyReci_LowThreads<16,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    8: EwaldGPU_EnergyReci_LowThreads<8,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    4: EwaldGPU_EnergyReci_LowThreads<4,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    2: EwaldGPU_EnergyReci_LowThreads<2,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
			case    1: EwaldGPU_EnergyReci_LowThreads<1,false><<<dimGrid6, dimBlock6, smemSize, *stream0>>>(m_dev_rho_hat,m_dev_infl_factor,m_dev_energy_reci,m_N,m_num_k,m_V,maxThreadsPerBlockEnergie,loops); break;
		}
	}

	//Copy the values back from the GPU to the CPU
	HANDLE_ERROR( cudaMemcpyAsync( m_energy_reci, m_dev_energy_reci,sizeof(real),cudaMemcpyDeviceToHost, *stream0 ) );
	m_energy_reci[0] = coulomb.prefactor * m_energy_reci[0];
}

//COMPUTE K's
void EwaldgpuForce::compute_num_k()
{
	int Index=0;//Arrayindex

	for (INT64 ix=0; ix<=m_num_kx; ix++)
	{
		for (INT64 iy=-m_num_ky; iy<=m_num_ky; iy++)
		{
			for (INT64 iz=-m_num_kz; iz<=m_num_kz; iz++)
			{
				if(ix*(2*m_num_ky+1)*(2*m_num_kz+1)+iy*(2*m_num_kz+1)+iz >= 0//Half m_k-space
				   and POW2(ix)*POW2((INT64)m_num_ky)*POW2((INT64)m_num_kz)
				   	 + POW2(iy)*POW2((INT64)m_num_kx)*POW2((INT64)m_num_kz)
				   	 + POW2(iz)*POW2((INT64)m_num_kx)*POW2((INT64)m_num_ky)
				   	 <=POW2((INT64)m_num_kx)*POW2((INT64)m_num_ky)*POW2((INT64)m_num_kz))//m_k-space ellipsoid
				{
					Index+=1;
				}
			}
		}
	}
	m_num_k=Index;
}
void EwaldgpuForce::compute_k_AND_influence_factor()
{
	real k_sqr;//Absolute square of m_k-vector
	int Index=0;//Arrayindex

	for (INT64 ix=0; ix<=m_num_kx; ix++)//Half m_k-space
	{
		for (INT64 iy=-m_num_ky; iy<=m_num_ky; iy++)
		{
			for (INT64 iz=-m_num_kz; iz<=m_num_kz; iz++)
			{
				if(ix*(2*m_num_ky+1)*(2*m_num_kz+1)+iy*(2*m_num_kz+1)+iz >= 0//Half m_k-space
				   and POW2(ix)*POW2((INT64)m_num_ky)*POW2((INT64)m_num_kz)
				   	 + POW2(iy)*POW2((INT64)m_num_kx)*POW2((INT64)m_num_kz)
				   	 + POW2(iz)*POW2((INT64)m_num_kx)*POW2((INT64)m_num_ky)
				   	 <=POW2((INT64)m_num_kx)*POW2((INT64)m_num_ky)*POW2((INT64)m_num_kz))//m_k-space ellipsoid
				{
					//m_k vectors
					m_k[Index] = 2*M_PI/m_box_l[0]*ix;
					m_k[Index+m_num_k] = 2*M_PI/m_box_l[1]*iy;
					m_k[Index+2*m_num_k] = 2*M_PI/m_box_l[2]*iz;
					//Influence factor
					k_sqr= POW2(m_k[Index]) + POW2(m_k[Index+m_num_k]) + POW2(m_k[Index+2*m_num_k]);
					m_infl_factor[Index] = 8*M_PI/m_V*expf(-k_sqr/(4*POW2(m_alpha)))/k_sqr;
					//Index
					Index+=1;
				}
			}
		}
	}	m_infl_factor[0] = 0;//Influence factor at m_k=(0,0,0)
}

//ELSE
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
void EwaldgpuForce::Transform_Vector3d_to_Array(SystemInterface &systemInterface)
{
	int Index = 0;
  for(SystemInterface::const_vec_iterator &it = systemInterface.rBegin(); it != systemInterface.rEnd(); ++it)
  {
			m_r_i[Index] =  (*it)[0];
      m_r_i[Index+m_N] = (*it)[1];
      m_r_i[Index+2*m_N] = (*it)[2];
      Index++;
  }
  Index = 0;
  for(SystemInterface::const_real_iterator &it = systemInterface.qBegin(); it != systemInterface.qEnd(); ++it)
  {
      m_q_i[Index] = *it;
      Index++;
  }

  Eigen::Vector3d box = systemInterface.box();
  m_box_l[0] = box(0);
  m_box_l[1] = box(1);
  m_box_l[2] = box(2);

}
void EwaldgpuForce::Transform_ForceArray_to_Vector3d()
{
	F.clear();
	Eigen::Vector3d Force;
  for(int i=0; i<m_N; i++)
  {
			Force[0]=coulomb.prefactor * (m_forces_reci[i]);
			Force[1]=coulomb.prefactor * (m_forces_reci[i+m_N]);
			Force[2]=coulomb.prefactor * (m_forces_reci[i+2*m_N]);
  	  F.push_back(Force);
  }
}
void EwaldgpuForce::Output()
{
	//Position
	int Index = 0;
	for(SystemInterface::const_vec_iterator &it = System->rBegin(); it != System->rEnd(); ++it)
	{
			printf("POS_IT: ID:%i  %f %f %f\n",Index,(*it)[0],(*it)[1],(*it)[2]);
			printf("POS_AR: ID:%i  %f %f %f\n",Index,m_r_i[Index],m_r_i[Index+m_N],m_r_i[Index+2*m_N]);
			Index++;
	}
	//Force
	for(int i=0; i<m_N; i++)
	{
			printf("FORCE_RECI_IT: ID:%i  %f %f %f\n",i,F[i][0],F[i][1],F[i][2]);
			printf("FORCE_RECI_AR: ID:%i  %f %f %f\n",i,m_forces_reci[i],m_forces_reci[i+m_N],m_forces_reci[i+2*m_N]);
	}
	//Energy
	printf("ENERGY_RECI_AR: %f\n",m_energy_reci[0]);
}

//REAL SPACE
void EwaldgpuForce::EwaldCPU_EnergySelf()
{
	m_energy_self[0]=coulomb.prefactor * (-m_alpha/sqrt(M_PI) * m_q_sqr);
}
void EwaldgpuForce::compute_q_sqare()
{
	m_q_sqr=0;
	for(int i=0;i<m_N;i++)
	{
		m_q_sqr += POW2(m_q_i[i]);
	}
}

//PARAMETERS
int ewaldgpu_set_params(double rcut, int num_kx, int num_ky, int num_kz, double alpha)
{
	IA_parameters *data = get_ia_param_safe(0, 0);
	data->max_cut = rcut;
	ewaldgpu_params.rcut = rcut;
	ewaldgpu_params.num_kx = num_kx;
	ewaldgpu_params.num_ky = num_ky;
	ewaldgpu_params.num_kz = num_kz;
	ewaldgpu_params.alpha = alpha;

	return 0;
}
int ewaldgpu_set_params_tune(double accuracy, double precision, int K_max, int time_calc_steps)
{
	ewaldgpu_params.accuracy = accuracy;
	ewaldgpu_params.precision = precision;
	ewaldgpu_params.K_max = K_max;
	ewaldgpu_params.time_calc_steps = time_calc_steps;

	return 0;
}

//TUNING r_cut, num_kx, num_ky, num_kz, alpha
int ewaldgpu_adaptive_tune(char **log)
{
	ewaldgpu_params.isTuned = false;
	int Kmax = ewaldgpu_params.K_max;
	double alpha_array[Kmax]; //All computed alpha in dependence of K
	double rcut_array[Kmax]; //All computed r_cut in dependence of all computed alpha
	double q_sqr = ewaldgpu_compute_q_sqare();
	char b[3*ES_INTEGER_SPACE + 3*ES_DOUBLE_SPACE + 128];

  if (skin == -1) {
    *log = strcat_alloc(*log, "ewaldgpu cannot be tuned, since the skin is not yet set");
    return ES_ERROR;
  }

	//Compute alpha for all reciprocal k-sphere radius K
	for(int K = 0; K < Kmax ;K++)
	{
		alpha_array[K] = ewaldgpu_tune_alpha(ewaldgpu_params.accuracy/sqrt(2), ewaldgpu_params.precision, K+1, box_l[0]*box_l[1]*box_l[2], q_sqr, n_total_particles);
		//printf("K:%i alpha:%f\n",K+1,alpha_array[K]);
	}
	//Compute r_cut for all computed alpha
	for(int K = 0; K < Kmax ;K++)
	{
		rcut_array[K] = ewaldgpu_tune_rcut(ewaldgpu_params.accuracy/sqrt(2), ewaldgpu_params.precision, alpha_array[K], box_l[0]*box_l[1]*box_l[2], q_sqr, n_total_particles);
		//printf("K:%i rcut:%f \n",K+1,rcut_array[K]);
	}
	//Test if accuracy was reached
	if(rcut_array[Kmax-1]<0)
  {
    return ES_ERROR;
  }

	/***********************************************************************************
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	PERFORMANCE TIME
	 ***********************************************************************************/

	//Test performance time for the diverent (K, rcut, alpha)
	double int_time_best = 1E30;
	int K_best;
	for(int K = 0; K < Kmax ;K++)
	{
		if(alpha_array[K]>0 and rcut_array[K]>0 and rcut_array[K]<(std::min(box_l[0],std::min(box_l[1],box_l[2])))/2.0-skin)
		{
			ewaldgpu_set_params(rcut_array[K], K+1, K+1, K+1, alpha_array[K]);
			mpi_bcast_coulomb_params();
			double int_time = time_force_calc(ewaldgpu_params.time_calc_steps);
			if(int_time<int_time_best)
			{
				int_time_best = int_time;
				K_best = K;
			}
			//printf("TIME K:%i int_time:%f\n",K+1,int_time);
		}
	}

	ewaldgpu_set_params(rcut_array[K_best], K_best+1, K_best+1, K_best+1, alpha_array[K_best]);
  ewaldgpu_params.isTuned = true;
	mpi_bcast_coulomb_params();

  /* Print Status */
  sprintf(b, "ewaldgpu tune parameters: Accuracy goal = %f\n", ewaldgpu_params.accuracy);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: Alpha = %f\n", ewaldgpu_params.alpha);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: r_cut = %f\n", ewaldgpu_params.rcut);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: num_kx = %i\n", ewaldgpu_params.num_kx);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: num_ky = %i\n", ewaldgpu_params.num_ky);
  *log = strcat_alloc(*log, b);
  sprintf(b, "ewaldgpu tune parameters: num_kz = %i\n", ewaldgpu_params.num_kz);
  *log = strcat_alloc(*log, b);

  return ES_OK;
}
double ewaldgpu_error_estimate_r(double q_sqr, int N, double r_cut, double V, double alpha, double accuracy)
{
	return sqrt(q_sqr/N) * 2 * (sqrt(q_sqr/(V*r_cut))) * exp(-alpha*alpha * r_cut*r_cut) - accuracy ;  // Kolafa-Perram, eq. 18
}
double ewaldgpu_error_estimate_k(double q_sqr, int N, int K, double V, double alpha, double accuracy)
{
	return sqrt(q_sqr/N) * alpha / (pow(V,1/3.0) * M_PI) * pow(8*q_sqr / K,0.5) * exp(-M_PI*M_PI * K*K / (alpha*alpha * pow(V,2/3.0))) - accuracy;  // Kolafa-Perram, eq. 32
}
double ewaldgpu_tune_alpha(double accuracy, double precision, int K, double V, double q_sqr, int N)
{
	double alpha_low=0.01;
	double alpha_high=100;
	double alpha_guess;
	double fkt_low;
	double fkt_high;
	double fkt_guess;

	// Find alpha with given K in k-space error estimate via bisection
	fkt_low = ewaldgpu_error_estimate_k(q_sqr, N, K, V, alpha_low, accuracy);
	fkt_high = ewaldgpu_error_estimate_k(q_sqr, N, K, V, alpha_high, accuracy);

  if (fkt_low*fkt_high > 0.0)
  {
    return -1; // Value unusable
  }

  do
  {
    alpha_guess = 0.5 *(alpha_low + alpha_high);
    fkt_guess = ewaldgpu_error_estimate_k(q_sqr, N, K, V, alpha_guess, accuracy);
    if (fkt_low*fkt_guess < 0.0) alpha_high = alpha_guess;
    else alpha_low = alpha_guess;

  } while (fabs(alpha_low-alpha_high) > precision);

  return 0.5 *(alpha_low + alpha_high);

}
double ewaldgpu_tune_rcut(double accuracy, double precision, double alpha, double V, double q_sqr, int N)
{
	double rcut_low=0.001;
	double rcut_high=0.5 * max(box_l[0],max(box_l[1],box_l[2]));
	double rcut_guess;
	double fkt_low;
	double fkt_high;
	double fkt_guess;

	// Find rcut with given K in k-space error estimate via bisection
	fkt_low = ewaldgpu_error_estimate_r(q_sqr,  N, rcut_low, V, alpha, accuracy);
	fkt_high = ewaldgpu_error_estimate_r(q_sqr, N, rcut_high, V, alpha, accuracy);

	if (fkt_low*fkt_high > 0.0)
  {
 		return -1; // Value unusable
  }

  do
  {
    rcut_guess = 0.5 *(rcut_low + rcut_high);
    fkt_guess = ewaldgpu_error_estimate_r(q_sqr, N, rcut_guess, V, alpha, accuracy);
    if (fkt_low*fkt_guess < 0.0) rcut_high = rcut_guess;
    else rcut_low = rcut_guess;

  } while (fabs(rcut_low-rcut_high) > precision);

  return 0.5 *(rcut_low + rcut_high);

}
int ewaldgpu_count_charged_particles()
{
  Cell *cell;
  Particle *part;
  int i,c,np;
  int sum_qpart =0;

  for (c = 0; c < local_cells.n; c++)
  {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++)
    {
      if( part[i].p.q != 0.0 )
      {
      	sum_qpart += 1.0;
      }
    }
  }
  sum_qpart    = (int)(sum_qpart+0.1);
  return (1999 + sum_qpart)/sum_qpart;
}

//KOLAFFA compute optimal alpha
double ewaldgpu_compute_E_error_estimate_r(double alpha, double rcut, double q_sqr, double box_l[3])
{
	//Compute the r space part of the force error estimate
  double std_deviation;
  std_deviation = q_sqr*pow(rcut/(2.0*box_l[0]*box_l[1]*box_l[2]),0.5) * exp(-pow(alpha,2)*pow(rcut,2)) / (pow(alpha*rcut,2));  // Kolafa-Perram, eq. 18

  return std_deviation;
}
double ewaldgpu_compute_E_error_estimate_k(double alpha, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3])
{
  //Compute the r space part of the force error estimate
  double std_deviation;
	std_deviation = q_sqr * alpha * pow(M_PI,-2.0) * pow(num_kx,-1.5) * exp(-(pow(M_PI*num_kx/(alpha*(box_l[0]+box_l[1]+box_l[2])/3.0),2))); //Kolafa Perram, eq. 32
  return std_deviation;
}
double ewaldgpu_E_estimate_error(double rcut, int num_kx, int num_ky, int num_kz, double alpha, double q_sqr, double box_l[3])
{
	//Compute the energy_reci error estimate
  return sqrt(pow(ewaldgpu_compute_E_error_estimate_r(alpha, rcut, q_sqr, box_l),2) + pow(ewaldgpu_compute_E_error_estimate_k(alpha, num_kx, num_ky, num_kz, q_sqr, box_l),2));
}
double ewaldgpu_compute_optimal_alpha(double rcut, int num_kx, int num_ky, int num_kz, double q_sqr, double box_l[3], double precision)
{
  //Use bisectional method to get optimal alpha value
  double alpha_low, f_low;
  double alpha_high, f_high;
  double alpha_guess, f_guess;
	int counter=0;
  alpha_low = 0.01;
  alpha_high = 10.0;

  f_low = ewaldgpu_compute_E_error_estimate_r(alpha_low, rcut, q_sqr, box_l) - ewaldgpu_compute_E_error_estimate_k(alpha_low, num_kx, num_ky, num_kz, q_sqr, box_l);
  f_high = ewaldgpu_compute_E_error_estimate_r(alpha_high, rcut, q_sqr, box_l) - ewaldgpu_compute_E_error_estimate_k(alpha_high, num_kx, num_ky, num_kz, q_sqr, box_l);

  if (f_low*f_high > 0.0)
  {
    fprintf(stderr, "Error: Could not init method to find optimal alpha!\n");
    exit(1);
  }

  do
  {
    alpha_guess = 0.5 *(alpha_low + alpha_high);
    f_guess = ewaldgpu_compute_E_error_estimate_r(alpha_guess, rcut, q_sqr, box_l) - ewaldgpu_compute_E_error_estimate_k(alpha_guess, num_kx, num_ky, num_kz, q_sqr, box_l);
    if (f_low*f_guess < 0.0) alpha_high = alpha_guess;
    else alpha_low = alpha_guess;
    counter++;
		if(counter>10000)
		{
			fprintf(stderr, "Find optimal alpha: Maximum number of loops exceeded: %i loops",counter);
			exit(1);
		}
  } while (fabs(alpha_low-alpha_high) > precision);

  return 0.5 *(alpha_low + alpha_high);
}
double ewaldgpu_compute_q_sqare()
{
	double q_sqr=0;
  Cell *cell;
  Particle *p;
  int i,c,np;

  for (c = 0; c < local_cells.n; c++)
  {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++)
    {
    	q_sqr += POW2(p[i].p.q);
    }
  }
  return q_sqr;
}
#endif
