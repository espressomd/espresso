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

#ifdef P3M
#ifdef CUDA

#include "electrostatics/p3m_gpu.hpp"
#include "electrostatics/p3m_gpu_cuda.cuh"

#include "actor/visitors.hpp"
#include "electrostatics/coulomb.hpp"

#include "ParticleRange.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"

#include "communication.hpp"

void CoulombP3MGPU::add_long_range_forces(ParticleRange const &) {
  if (this_node == 0) {
    p3m_gpu_add_farfield_force(prefactor);
  }
}

void CoulombP3MGPU::init() {
  if (has_actor_of_type<ElectrostaticLayerCorrection>(electrostatics_actor)) {
    init_cpu_kernels();
  }
  p3m_gpu_init(p3m.params.cao, p3m.params.mesh.data(), p3m.params.alpha);
}

void CoulombP3MGPU::init_cpu_kernels() { CoulombP3M::init(); }

void CoulombP3MGPU::request_gpu() const {
  auto &gpu_particle_data = System::get_system().gpu;
  gpu_particle_data.enable_property(GpuParticleData::prop::force);
  gpu_particle_data.enable_property(GpuParticleData::prop::q);
  gpu_particle_data.enable_property(GpuParticleData::prop::pos);
}

#endif // CUDA
#endif // P3M
