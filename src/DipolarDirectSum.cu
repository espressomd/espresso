#include "config.hpp"

#ifdef DIPOLAR_DIRECT_SUM
#include "DipolarDirectSum.hpp"
#include "grid.hpp"

#include "cuda_utils.hpp"
#include <stdio.h>
#include "dipole_cuda.hpp"



DipolarDirectSum *dipolarDirectSum = 0;

__global__ void DipolarDirectSum_kernel(float k,
				     int n, float *pos, float* dip, float *f, float* torque, float box_l[3], int periodic[3]) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  float pf;

  if(i >= n)
    return;
  
  for (int j=0;j<n;j++)
  {
   pf = (i!=j) *k;
   dipole_ia(pf,pos+3*i,pos+3*j,dip+3*i,dip+3*j,f+3*i,torque+3*i,1,box_l,periodic);
 }
}


void DipolarDirectSum_kernel_wrapper(float k, int n, float *pos, float *dip, float* f, float* torque) {
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


  float bl[3];
  int per[3];
  for (int i=0;i<3;i++)
  {
   bl[i]=box_l[i];
  #ifdef PARTIAL_PERIODIC
   per[i]=PERIODIC(i);
  #else
   per[i]=1;
  #endif
 }

  KERNELCALL(DipolarDirectSum_kernel,grid,block,(k, n, pos, dip,f,torque,bl,per));
}
#endif
