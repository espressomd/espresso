
#include "integrate_sd_cuda.cuh"
#include "integrate_sd_cuda_device.cuh"

/* *************************************************************************************************************** *
 * ********************************************    DEVICE-Functions   ******************************************** *
 * *************************************************************************************************************** */

__device__ double atomicAdd(double * address, double inc){

  ull *addressUll = (ull*) address;
  ull oldValue=*addressUll;
  ull assumedValue;
  assert(!isnan(inc));
  do {
    assumedValue=oldValue;
    ull newValue = __double_as_longlong (__longlong_as_double(assumedValue)+inc);
    oldValue = atomicCAS(addressUll,assumedValue,newValue);
  }
  while (oldValue != assumedValue);
  return __longlong_as_double(oldValue);
}

__device__ __inline__ int fold_back_up(int value, int max){
  while (value < 0){
    value+=max;
  }
  while (value >= max){
    value-=max;
  }
  return value;
}

__device__ int reduce_max(int * shared_cache, int value){
  shared_cache[threadIdx.x]=value;
  shared_cache[threadIdx.x+blockDim.x]=0;
  for (int t=(blockDim.x+1)/2;t>1;t=(t+1)/2){
    if (threadIdx.x < t){
      shared_cache[threadIdx.x]=shared_cache[threadIdx.x]>shared_cache[threadIdx.x+t]?shared_cache[threadIdx.x]:shared_cache[threadIdx.x+t];
    }
    __syncthreads();
  }
  if (threadIdx.x==0){
    shared_cache[0]=shared_cache[0]>shared_cache[1]?shared_cache[0]:shared_cache[1];
  }
  return shared_cache[0];
}
__device__ void reduce_sum(real * shared_cache){
#ifndef numThreadsPerBlock_is_power_of_two
#error numThreadsPerBlock has to be a power of two for effective reduction
#endif
  //shared_cache[threadIdx.x]=value;
  for (int t=(blockDim.x)/2;t>0;t=t/2){
    if (threadIdx.x < t){
      shared_cache[threadIdx.x]=shared_cache[threadIdx.x]+shared_cache[threadIdx.x+t];
    }
    __syncthreads();
  }
}

__device__ __inline__ real read_without_caching( const real * addr){
   real tmp;
#ifdef SD_USE_FLOAT
   asm("ld.global.cs.f32 %0,[%1];\n"
       : "=f"(tmp) : "l"(addr) : );
#else
   asm("ld.global.cs.f64 %0,[%1];\n"
       : "=d"(tmp) : "l"(addr) : );
#endif
   return tmp;
}
