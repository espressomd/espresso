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
#include "EspressoSystemInterface.hpp"
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "grid.hpp"

/* Initialize instance pointer */
EspressoSystemInterface *EspressoSystemInterface::m_instance = nullptr;

/********************************************************************************************/

void EspressoSystemInterface::gatherParticles() {
// get particles from other nodes
#ifdef CUDA
  if (m_gpu) {
    if (gpu_get_global_particle_vars_pointer_host()->communication_enabled) {
      copy_part_data_to_gpu(cell_structure.local_particles());
      reallocDeviceMemory(gpu_get_particle_pointer().size());
      if (m_splitParticleStructGpu && (this_node == 0))
        split_particle_struct();
    }
  }
#endif
}

void EspressoSystemInterface::init() { gatherParticles(); }

void EspressoSystemInterface::update() { gatherParticles(); }

Utils::Vector3d EspressoSystemInterface::box() const {
  return box_geo.length();
}
