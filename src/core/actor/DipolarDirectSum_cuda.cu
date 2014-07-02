
#include "config.hpp"

#ifdef DIPOLAR_DIRECT_SUM
#include "DipolarDirectSum.hpp"

#include "cuda_utils.hpp"
#include <stdio.h>
#include "dipole_cuda.hpp"




__global__ void DipolarDirectSum_kernel(float k,
				     int n, float *pos, float* dip, float *f, float* torque, float box_l[3], int periodic[3]) {

  printf("kkkk-box_l: %f %f %f\n",box_l[0],box_l[1],box_l[2]);
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  float pf;

printf("%d: %f %f %f\n",i,*(f+3*i),*(f+3*i+1),*(f+3*i+2));
  if(i >= n)
    return;
  
  for (int j=0;j<n;j++)
  {

   printf("%d %d\n",i,j);
   if (i!=j)
     dipole_ia(1,pos+3*i,pos+3*j,dip+3*i,dip+3*j,f+3*i,torque+3*i,1,box_l,periodic);

 }
printf("%d: %f %f %f\n",i,*(f+3*i),*(f+3*i+1),*(f+3*i+2));
}


void DipolarDirectSum_kernel_wrapper(float k, int n, float *pos, float *dip, float* f, float* torque, float box_l[3],int periodic[3]) {
  dim3 grid(1,1,1);
  dim3 block(1,1,1);

  if(n == 0)
    return;

  if(n <= 512) {
    grid.x = 1;
    block.x = n;
  } else {
    grid.x = n/512 + 1;
    block.x = 512;
  }

  float* box_l_gpu;
  int* periodic_gpu;
  cuda_safe_mem(cudaMalloc((void**)&box_l_gpu,3*sizeof(float)));
  cuda_safe_mem(cudaMalloc((void**)&periodic_gpu,3*sizeof(int)));
  cuda_safe_mem(cudaMemcpy(box_l_gpu,box_l,3*sizeof(float),cudaMemcpyHostToDevice));
  cuda_safe_mem(cudaMemcpy(periodic_gpu,periodic,3*sizeof(int),cudaMemcpyHostToDevice));



  printf("box_l: %f %f %f\n",box_l[0],box_l[1],box_l[2]);
  KERNELCALL(DipolarDirectSum_kernel,grid,block,(k, n, pos, dip,f,torque,box_l_gpu, periodic_gpu));
  cudaFree(box_l_gpu);
  cudaFree(periodic_gpu);

}


#endif
