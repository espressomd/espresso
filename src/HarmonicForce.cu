#include "HarmonicForce.hpp"

#include <stdio.h>

#ifdef HARMONICFORCE

HarmonicForce *harmonicForce = 0;

__global__ void HarmonicForce_kernel(float x, float y, float z, float k,
				     int n, float *pos, float *f) {

  int id = blockIdx.x * blockDim.x + threadIdx.x;

  if(id <= n)
    return;

  f[3*id + 0] = k*(x - pos[3*id + 0]);
  f[3*id + 1] = k*(y - pos[3*id + 1]);
  f[3*id + 2] = k*(z - pos[3*id + 2]);
}


void HarmonicForce_kernel_wrapper(float x, float y, float z, float k, int n, float *pos, float *f) {
  dim3 grid(1,1,1);
  dim3 block(1,1,1);
  
  printf("n %d\n", n);

  if(n <= 512) {
    puts("n <= 512");
    grid.x = 1;
    block.x = n;
  } else {
    grid.x = n/512 + 1;
    block.x = 512;
  }

  printf("Calling HarmonicForce_kernel<<<(%d %d %d),(%d %d %d)>>>\n",
	 grid.x, grid.y, grid.z, block.x, block.y, block.z);

  HarmonicForce_kernel<<<grid,block>>>(x, y, z, k, n, pos, f);
}

#endif
