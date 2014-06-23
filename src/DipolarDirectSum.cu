#include "DipolarDirectSum.hpp"

#include "cuda_utils.hpp"
#include <stdio.h>

#ifdef DIPOLAR_DIRECT_SUM

DipolarDirectSum *dipolarDirectSum = 0;

__global__ void dipolarDirectSum_kernel(float k,
				     int n, float *pos, float* dip, float *f, float* torque) {

  int id = blockIdx.x * blockDim.x + threadIdx.x;
  float pf;

  if(id >= n)
    return;
  
  for (j=0;j<n;j++)
  {
   pf = (i!=j) *k;
   dipole_ia(pf,pos[id],pos[j],dip[id],dip[j],f[id],torque[id],1);
}


void DipolarDirectSum_kernel_wrapper(float x, float y, float z, float k, int n, float *pos, float *f, float *dip) {
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

  KERNELCALL(DipolarDirectSum_kernel,grid,block,(x, y, z, k, n, pos, f,dip))
}
#endif
