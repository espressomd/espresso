/*
  Copyright (C) 2014-2018 The ESPResSo project

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

#define _P3M_GPU_FLOAT
//#define _P3M_GPU_REAL_DOUBLE

#ifdef _P3M_GPU_FLOAT
#define REAL_TYPE float
#define FFT_TYPE_COMPLEX hipfftComplex
#define FFT_FORW_FFT hipfftExecR2C
#define FFT_BACK_FFT hipfftExecC2R
#define FFT_PLAN_FORW_FLAG HIPFFT_R2C
#define FFT_PLAN_BACK_FLAG HIPFFT_C2R
#endif

#ifdef _P3M_GPU_REAL_DOUBLE
#define REAL_TYPE double
#define FFT_TYPE_COMPLEX hipfftDoubleComplex
#define FFT_FORW_FFT hipfftExecD2Z
#define FFT_BACK_FFT hipfftExecZ2D
#define FFT_PLAN_FORW_FLAG HIPFFT_D2Z
#define FFT_PLAN_BACK_FLAG HIPFFT_Z2D
#endif

void p3m_gpu_init(int cao, int mesh[3], double alpha);
void p3m_gpu_add_farfield_force();

#endif /* _P3M_GPU_H */
