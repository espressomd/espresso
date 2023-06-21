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

#include "communication.hpp"
#include "grid.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"

static void get_simulation_box(float *box, int *per) {
  for (int i = 0; i < 3; i++) {
    box[i] = static_cast<float>(box_geo.length()[i]);
    per[i] = box_geo.periodic(i);
  }
}

DipolarDirectSumGpu::DipolarDirectSumGpu(double prefactor)
    : prefactor{prefactor} {
  auto &gpu_particle_data = System::get_system().gpu;
  gpu_particle_data.enable_property(GpuParticleData::prop::force);
  gpu_particle_data.enable_property(GpuParticleData::prop::torque);
  gpu_particle_data.enable_property(GpuParticleData::prop::pos);
  gpu_particle_data.enable_property(GpuParticleData::prop::dip);
}

void DipolarDirectSumGpu::add_long_range_forces() const {
  auto &gpu = System::get_system().gpu;
  gpu.update();
  if (this_node != 0) {
    return;
  }
  float box[3];
  int periodicity[3];
  get_simulation_box(box, periodicity);
  auto const npart = static_cast<unsigned>(gpu.n_particles());
  auto const forces_device = gpu.get_particle_forces_device();
  auto const torques_device = gpu.get_particle_torques_device();
  auto const positions_device = gpu.get_particle_positions_device();
  auto const dipoles_device = gpu.get_particle_dipoles_device();
  DipolarDirectSum_kernel_wrapper_force(
      static_cast<float>(prefactor), npart, positions_device, dipoles_device,
      forces_device, torques_device, box, periodicity);
}

void DipolarDirectSumGpu::long_range_energy() const {
  auto &gpu = System::get_system().gpu;
  gpu.update();
  if (this_node != 0) {
    return;
  }
  float box[3];
  int periodicity[3];
  get_simulation_box(box, periodicity);
  auto const npart = static_cast<unsigned>(gpu.n_particles());
  auto const energy_device = &(gpu.get_energy_device()->dipolar);
  auto const positions_device = gpu.get_particle_positions_device();
  auto const dipoles_device = gpu.get_particle_dipoles_device();
  DipolarDirectSum_kernel_wrapper_energy(static_cast<float>(prefactor), npart,
                                         positions_device, dipoles_device, box,
                                         periodicity, energy_device);
}

#endif
