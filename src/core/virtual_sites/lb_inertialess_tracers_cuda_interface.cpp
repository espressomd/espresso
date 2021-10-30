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

// This is an internal file of the IMMERSED BOUNDARY implementation
// It should not be included by any main ESPResSo routines
// Functions to be exported for ESPResSo are in ibm_main.hpp

#include "config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

#include "Particle.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "serialization/ibm_cuda_particle_velocities_input.hpp"
#include "virtual_sites/lb_inertialess_tracers_cuda_interface.hpp"

#include <utils/mpi/gather_buffer.hpp>
#include <utils/mpi/scatter_buffer.hpp>

#include <vector>

// Variables for communication
std::vector<IBM_CUDA_ParticleDataInput> IBM_ParticleDataInput_host = {};
std::vector<IBM_CUDA_ParticleDataOutput> IBM_ParticleDataOutput_host = {};

static void pack_particles(ParticleRange const &particles,
                           std::vector<IBM_CUDA_ParticleDataInput> &buffer) {

  int i = 0;
  for (auto const &part : particles) {
    auto const pos = folded_position(part.r.p, box_geo);

    buffer[i].pos[0] = static_cast<float>(pos[0]);
    buffer[i].pos[1] = static_cast<float>(pos[1]);
    buffer[i].pos[2] = static_cast<float>(pos[2]);

    buffer[i].f[0] = static_cast<float>(part.f.f[0]);
    buffer[i].f[1] = static_cast<float>(part.f.f[1]);
    buffer[i].f[2] = static_cast<float>(part.f.f[2]);

    buffer[i].is_virtual = part.p.is_virtual;

    i++;
  }
}

/** Gather particle positions on the head node in order to communicate them
 *  to GPU. We transfer all particles (real and virtual), but actually we would
 *  only need the virtual ones. Room for improvement...
 *  Analogous to @ref cuda_mpi_get_particles.
 */
void IBM_cuda_mpi_get_particles(ParticleRange const &particles) {
  auto const n_part = particles.size();

  if (this_node > 0) {
    static std::vector<IBM_CUDA_ParticleDataInput> buffer;
    buffer.resize(n_part);
    /* pack local parts into buffer */
    pack_particles(particles, buffer);

    Utils::Mpi::gather_buffer(buffer, comm_cart);
  } else {
    /* Pack own particles */
    pack_particles(particles, IBM_ParticleDataInput_host);

    Utils::Mpi::gather_buffer(IBM_ParticleDataInput_host, comm_cart);
  }
}

static void set_velocities(ParticleRange const &particles,
                           std::vector<IBM_CUDA_ParticleDataOutput> &buffer) {
  int i = 0;
  for (auto &part : particles) {
    if (part.p.is_virtual) {
      for (int j = 0; j < 3; j++)
        part.m.v[j] = static_cast<double>(buffer[i].v[j]);
    }
    i++;
  }
}

/** Particle velocities have been communicated from GPU, now transmit to all
 *  nodes. Analogous to @ref cuda_mpi_send_forces.
 */
void IBM_cuda_mpi_send_velocities(ParticleRange const &particles) {
  auto const n_part = particles.size();

  if (this_node > 0) {
    static std::vector<IBM_CUDA_ParticleDataOutput> buffer;
    /* Alloc buffer */
    buffer.resize(n_part);

    Utils::Mpi::scatter_buffer(buffer.data(), n_part, comm_cart);

    set_velocities(particles, buffer);
  } else {
    /* Scatter forces */
    Utils::Mpi::scatter_buffer(IBM_ParticleDataOutput_host.data(), n_part,
                               comm_cart);

    set_velocities(particles, IBM_ParticleDataOutput_host);
  }
}

#endif
