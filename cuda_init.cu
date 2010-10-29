#include <cuda.h>

extern "C" {

#include "utils.h"

}

#include "cuda_init.h"

void gpu_init()
{
  int deviceCount, dev, found;
  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess) {
    fprintf(stderr, "%d: cannot start CUDA.\n", this_node);
    errexit();
  }
  if (deviceCount == 0) {
    fprintf(stderr, "%d: There is no CUDA device.\n", this_node);
    errexit();
  }
  found = 0;
  for (dev = 0; dev < deviceCount; ++dev) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    if (deviceProp.major > 1 || (deviceProp.major == 1 && deviceProp.minor >= 1)) {
      fprintf(stderr, "%d:using CUDA device %s.\n", this_node, deviceProp.name);
      found = 1; break;
    }
  }
  if (!found) {
    fprintf(stderr, "%d: There is no device with compute capability >= 1.1.\n", this_node);
    errexit();
  }
}

