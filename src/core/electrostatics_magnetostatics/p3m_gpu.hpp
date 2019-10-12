/*
 * Copyright (C) 2014-2019 The ESPResSo project
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
#ifndef _P3M_GPU_H
#define _P3M_GPU_H
/** \file
 *  P3M electrostatics on GPU.
 *
 *  Implementation in p3m_gpu_cuda.cu.
 */

void p3m_gpu_init(int cao, const int mesh[3], double alpha);
void p3m_gpu_add_farfield_force();

#endif /* _P3M_GPU_H */
