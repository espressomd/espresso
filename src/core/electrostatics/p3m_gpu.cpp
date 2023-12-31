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
#include "PropagationMode.hpp"
#include "integrators/Propagation.hpp"
#include "npt.hpp"
#include "system/GpuParticleData.hpp"
#include "system/System.hpp"

#include "communication.hpp"

#include <cassert>
#include <limits>

static auto get_n_part_safe(GpuParticleData const &gpu) {
  auto const n_part = gpu.n_particles();
#ifndef NDEBUG
  auto constexpr n_part_max = std::numeric_limits<unsigned int>::max();
  assert(n_part < static_cast<std::size_t>(n_part_max));
#endif
  return static_cast<unsigned int>(n_part);
}

void CoulombP3MGPU::add_long_range_forces(ParticleRange const &particles) {
#ifdef NPT
  if (get_system().propagation->integ_switch == INTEG_METHOD_NPT_ISO) {
    auto const energy = long_range_energy(particles);
    npt_add_virial_contribution(energy);
  }
#else
  static_cast<void>(particles);
#endif
  if (this_node == 0) {
    auto &gpu = get_system().gpu;
    p3m_gpu_add_farfield_force(*m_gpu_data, gpu, prefactor,
                               get_n_part_safe(gpu));
  }
}

void CoulombP3MGPU::init() {
  auto &system = get_system();
  if (has_actor_of_type<ElectrostaticLayerCorrection>(
          system.coulomb.impl->solver)) {
    init_cpu_kernels();
  }
  p3m_gpu_init(m_gpu_data, p3m.params.cao, p3m.params.mesh.data(),
               p3m.params.alpha, system.box_geo->length(),
               get_n_part_safe(system.gpu));
}

void CoulombP3MGPU::init_cpu_kernels() { CoulombP3M::init(); }

void CoulombP3MGPU::request_gpu() const {
  auto &gpu_particle_data = get_system().gpu;
  gpu_particle_data.enable_property(GpuParticleData::prop::force);
  gpu_particle_data.enable_property(GpuParticleData::prop::q);
  gpu_particle_data.enable_property(GpuParticleData::prop::pos);
}

#endif // CUDA
#endif // P3M
