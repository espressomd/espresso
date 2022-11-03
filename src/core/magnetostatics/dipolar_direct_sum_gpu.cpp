/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef DIPOLAR_DIRECT_SUM

#include "magnetostatics/dipolar_direct_sum_gpu.hpp"
#include "magnetostatics/dipolar_direct_sum_gpu_cuda.cuh"

#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "grid.hpp"

static void get_simulation_box(float *box, int *per) {
  for (int i = 0; i < 3; i++) {
    box[i] = static_cast<float>(box_geo.length()[i]);
    per[i] = box_geo.periodic(i);
  }
}

DipolarDirectSumGpu::DipolarDirectSumGpu(double prefactor)
    : prefactor{prefactor} {
  auto &system = EspressoSystemInterface::Instance();
  system.requestFGpu();
  system.requestTorqueGpu();
  system.requestRGpu();
  system.requestDipGpu();
}

void DipolarDirectSumGpu::add_long_range_forces() const {
  auto &system = EspressoSystemInterface::Instance();
  system.update();
  if (this_node != 0) {
    return;
  }
  float box[3];
  int periodicity[3];
  get_simulation_box(box, periodicity);
  DipolarDirectSum_kernel_wrapper_force(
      static_cast<float>(prefactor), static_cast<unsigned>(system.npart_gpu()),
      system.rGpuBegin(), system.dipGpuBegin(), system.fGpuBegin(),
      system.torqueGpuBegin(), box, periodicity);
}

void DipolarDirectSumGpu::long_range_energy() const {
  auto &system = EspressoSystemInterface::Instance();
  system.update();
  if (this_node != 0) {
    return;
  }
  float box[3];
  int periodicity[3];
  get_simulation_box(box, periodicity);
  auto energy = &(reinterpret_cast<CUDA_energy *>(system.eGpu())->dipolar);
  DipolarDirectSum_kernel_wrapper_energy(
      static_cast<float>(prefactor), static_cast<unsigned>(system.npart_gpu()),
      system.rGpuBegin(), system.dipGpuBegin(), box, periodicity, energy);
}

#endif
