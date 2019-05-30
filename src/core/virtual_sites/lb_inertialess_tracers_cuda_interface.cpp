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

// *******
// This is an internal file of the IMMERSED BOUNDARY implementation
// It should not be included by any main Espresso routines
// Functions to be exported for Espresso are in ibm_main.hpp

#include "config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

#include "communication.hpp"
#include "debug.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "serialization/ibm_cuda_particle_velocities_input.hpp"
#include "virtual_sites/lb_inertialess_tracers_cuda_interface.hpp"

#include <utils/mpi/gather_buffer.hpp>
#include <utils/mpi/scatter_buffer.hpp>

// Variables for communication
IBM_CUDA_ParticleDataInput *IBM_ParticleDataInput_host = nullptr;
IBM_CUDA_ParticleDataOutput *IBM_ParticleDataOutput_host = nullptr;

/*****************
   IBM_cuda_mpi_get_particles
Gather particle positions on the master node in order to communicate them to GPU
We transfer all particles (real and virtual), but actually we would only need
the virtual ones
Room for improvement...
 *****************/

namespace {
void pack_particles(ParticleRange particles,
                    IBM_CUDA_ParticleDataInput *buffer) {
  int dummy[3] = {0, 0, 0};

  int i = 0;
  for (auto const &part : particles) {
    Utils::Vector3d pos = folded_position(part);

    buffer[i].pos[0] = (float)pos[0];
    buffer[i].pos[1] = (float)pos[1];
    buffer[i].pos[2] = (float)pos[2];

    buffer[i].f[0] = (float)part.f.f[0];
    buffer[i].f[1] = (float)part.f.f[1];
    buffer[i].f[2] = (float)part.f.f[2];

    buffer[i].is_virtual = part.p.is_virtual;

    i++;
  }
}
} // namespace

// Analogous to the usual cuda_mpi_get_particles function
void IBM_cuda_mpi_get_particles(ParticleRange particles) {
  auto const n_part = particles.size();

  if (this_node > 0) {
    COMM_TRACE(fprintf(stderr, "%d: get_particles_slave, %d particles\n",
                       this_node, n_part));
    static std::vector<IBM_CUDA_ParticleDataInput> buffer;
    buffer.resize(n_part);
    /* pack local parts into buffer */
    pack_particles(particles, buffer.data());

    Utils::Mpi::gather_buffer(buffer.data(), buffer.size(), comm_cart);
  } else {
    /* Pack own particles */
    pack_particles(particles, IBM_ParticleDataInput_host);

    Utils::Mpi::gather_buffer(IBM_ParticleDataInput_host, n_part, comm_cart);
  }

  COMM_TRACE(fprintf(stderr, "%d: finished get\n", this_node));
}

/*****************
   IBM_cuda_mpi_send_velocities
Particle velocities have been communicated from GPU, now transmit to all nodes
 ******************/
// Analogous to cuda_mpi_send_forces

namespace {
void set_velocities(ParticleRange particles,
                    IBM_CUDA_ParticleDataOutput *buffer) {
  int i = 0;
  for (auto &part : particles) {
    if (part.p.is_virtual)
      for (int j = 0; j < 3; j++)
        part.m.v[j] = buffer[i].v[j];

    i++;
  }
}
} // namespace

void IBM_cuda_mpi_send_velocities(ParticleRange particles) {
  auto const n_part = particles.size();

  if (this_node > 0) {
    static std::vector<IBM_CUDA_ParticleDataOutput> buffer;
    /* Alloc buffer */
    buffer.resize(n_part);

    Utils::Mpi::scatter_buffer(buffer.data(), n_part, comm_cart);

    set_velocities(particles, buffer.data());
  } else {
    /* Scatter forces to slaves */
    Utils::Mpi::scatter_buffer(IBM_ParticleDataOutput_host, n_part, comm_cart);

    set_velocities(particles, IBM_ParticleDataOutput_host);
  }
}

#endif
