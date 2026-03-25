/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include <config/config.hpp>

#if defined(ESPRESSO_DIPOLES) and defined(ESPRESSO_CUDA)

void DipolarDirectSum_kernel_wrapper_energy(float k, unsigned int n,
                                            float const *pos, float const *dip,
                                            float box_l[3], int periodic[3],
                                            float *E);
void DipolarDirectSum_kernel_wrapper_force(float k, unsigned int n,
                                           float const *pos, float const *dip,
                                           float *dip_fld, float *f,
                                           float *torque, float box_l[3],
                                           int periodic[3], int n_replicas);

#endif // defined(ESPRESSO_DIPOLES) and defined(ESPRESSO_CUDA)
