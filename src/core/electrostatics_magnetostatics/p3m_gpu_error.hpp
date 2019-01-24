/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef P3M_ERROR_GPU_HPP
#define P3M_ERROR_GPU_HPP
/** @file
 *  P3M electrostatics on GPU.
 *
 *  Implementation in p3m_gpu_error_cuda.cu.
 */

#include "config.hpp"

#ifdef CUDA
double p3m_k_space_error_gpu(double prefactor, const int *mesh, int cao,
                             int npart, double sum_q2, double alpha_L,
                             double *box);
#endif
#endif
