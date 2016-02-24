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

#ifndef __INTEGRATE_SD_CUDA_DEVICE_HPP
#define __INTEGRATE_SD_CUDA_DEVICE_HPP

#include "integrate_sd_cuda.hpp"

#ifdef SD

typedef unsigned long long ull;
// atomicAdd implementation for double
__device__ double atomicAdd(double * address, double inc);

__device__ int fold_back_up(int value, int max);

__device__ int reduce_max(int * shared_cache, int value);

__device__ __inline__ real read_without_caching( const real * addr);

__device__ void reduce_sum(real * shared_cache);

#endif /* SD */

#endif
