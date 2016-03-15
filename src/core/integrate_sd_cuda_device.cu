/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "integrate_sd_cuda.hpp"
#include "integrate_sd_cuda_device.hpp"

#ifdef SD

/* *************************************************************************************************************** *
 * ********************************************    DEVICE-Functions   ******************************************** *
 * *************************************************************************************************************** */

/// atomic add for doubles
/// address  : pointer to which the value should be added
/// inc      : value which should be added
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

/// modulo function for integers as implemented in python
/// returns value % max 
/// which is always positive
__device__ __inline__ int fold_back_up(int value, int max){
  while (value < 0){
    value+=max;
  }
  while (value >= max){
    value-=max;
  }
  return value;
}

/// reduction function returning maximum of all value.
/// has to be called by all threads
/// shared_cache should be a pointer to shared memory of at least size blockDim.x
/// blockDim.x has to be an even number
__device__ int reduce_max(int * shared_cache){
  for (int t=(blockDim.x+1)/2;t>1;t=(t+1)/2){
    if (threadIdx.x < t){
      if (shared_cache[threadIdx.x]<shared_cache[threadIdx.x+t]){
	shared_cache[threadIdx.x]=shared_cache[threadIdx.x+t];
      }
    }
    __syncthreads();
  }
  if (threadIdx.x==0){
    shared_cache[0]=shared_cache[0]>shared_cache[1]?shared_cache[0]:shared_cache[1];
  }
  return shared_cache[0];
}


__device__ int reduce_max(int * shared_cache, int value){
  shared_cache[threadIdx.x]=value;
  /*if (threadIdx.x == 0){
    if ((((blockDim.x+1)/2)*2) > blockDim.x){
      shared_cache[blockDim.x+1]=0;
    }
  }*/
  //shared_cache[threadIdx.x+blockDim.x]=0;
  __syncthreads();
  return reduce_max(shared_cache);
}



/// reduction function returning sum of all values in shared_cache[0:blockDim.x-1]
/// has to be called by all threads
/// shared_cache should be a pointer to shared memory of at least size blockDim.x
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

/// function to avoid caching if reading from global memory
/// addr     : address from which the value should be read
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

#endif
