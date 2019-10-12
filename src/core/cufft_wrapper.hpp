/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CUFFT_WRAPPER_HPP
#define CUFFT_WRAPPER_HPP

#if defined(__HIPCC__) // AMD or Nvidia via HIP

#include <hipfft.h>

#define cufftComplex hipfftComplex
#define cufftDestroy hipfftDestroy
#define cufftDoubleComplex hipfftDoubleComplex
#define cufftExecC2R hipfftExecC2R
#define cufftExecD2Z hipfftExecD2Z
#define cufftExecR2C hipfftExecR2C
#define cufftExecZ2D hipfftExecZ2D
#define cufftHandle hipfftHandle
#define cufftPlan3d hipfftPlan3d
#define cufftReal hipfftReal
#define cufftSetStream hipfftSetStream
#define CUFFT_R2C HIPFFT_R2C
#define CUFFT_C2R HIPFFT_C2R
#define CUFFT_D2Z HIPFFT_D2Z
#define CUFFT_Z2D HIPFFT_Z2D
#define CUFFT_SUCCESS HIPFFT_SUCCESS

#else // Nvidia via CUDA

#include <cufft.h>

#endif

#endif // CUFFT_WRAPPER_HPP
