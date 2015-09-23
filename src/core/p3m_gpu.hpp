/*
  Copyright (C) 2014 The ESPResSo project
  
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
#ifndef _P3M_GPU_H
#define _P3M_GPU_H

//NOTE :if one wants to use doubles it requires cuda compute capability 1.3
#define _P3M_GPU_FLOAT
//#define _P3M_GPU_REAL_DOUBLE

#ifdef _P3M_GPU_FLOAT
#define REAL_TYPE float
#define CUFFT_TYPE_COMPLEX cufftComplex
#define CUFFT_FFT cufftExecC2C
#define CUFFT_PLAN_FLAG CUFFT_C2C
#endif

#ifdef _P3M_GPU_REAL_DOUBLE
#define REAL_TYPE double
#define CUFFT_TYPE_COMPLEX cufftDoubleComplex
#define CUFFT_FFT cufftExecZ2Z
#define CUFFT_PLAN_FLAG CUFFT_Z2Z
#endif

void p3m_gpu_init(int cao, int mesh[3], double alpha, double box[3]);
void p3m_gpu_add_farfield_force();

#endif /* _P3M_GPU_H */

