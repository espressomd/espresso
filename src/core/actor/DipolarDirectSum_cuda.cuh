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

#ifndef DIPOLARDIRECTSUM_CUH_
#define DIPOLARDIRECTSUM_CUH_

#include "config.hpp"

#ifdef DIPOLAR_DIRECT_SUM

void DipolarDirectSum_kernel_wrapper_energy(float k, unsigned int n, float *pos,
                                            float *dip, float box_l[3],
                                            int periodic[3], float *E);
void DipolarDirectSum_kernel_wrapper_force(float k, unsigned int n, float *pos,
                                           float *dip, float *f, float *torque,
                                           float box_l[3], int periodic[3]);

#endif // DIPOLAR_DIRECT_SUM

#endif /* DIPOLARDIRECTSUM_CUH_ */
