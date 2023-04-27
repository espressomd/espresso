/*
 * Copyright (C) 2014-2022 The ESPResSo project
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

#include "cuda_interface.hpp"

#ifdef CUDA

#include "EspressoSystemInterface.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "grid.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "serialization/CUDA_particle_data.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>
#include <utils/mpi/scatter_buffer.hpp>

#include <vector>

static void pack_particles(ParticleRange particles,
                           CUDA_particle_data *buffer) {
  using Utils::Vector3f;

  int i = 0;
  for (auto const &part : particles) {
    buffer[i].p = static_cast<Vector3f>(folded_position(part.pos(), box_geo));

    buffer[i].identity = part.id();
    buffer[i].v = static_cast<Vector3f>(part.v());
#ifdef VIRTUAL_SITES
    buffer[i].is_virtual = part.is_virtual();
#endif

#ifdef DIPOLES
    buffer[i].dip = static_cast<Vector3f>(part.calc_dip());
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
    buffer[i].mu_E = static_cast<Vector3f>(part.mu_E());
#endif

#ifdef ELECTROSTATICS
    buffer[i].q = static_cast<float>(part.q());
#endif

#ifdef MASS
    buffer[i].mass = static_cast<float>(part.mass());
#endif

#ifdef ROTATION
    buffer[i].director = static_cast<Vector3f>(part.calc_director());
#endif

#ifdef ENGINE
    buffer[i].swim.v_swim = static_cast<float>(part.swimming().v_swim);
    buffer[i].swim.f_swim = static_cast<float>(part.swimming().f_swim);
    buffer[i].swim.director = buffer[i].director;

    buffer[i].swim.push_pull = part.swimming().push_pull;
    buffer[i].swim.dipole_length =
        static_cast<float>(part.swimming().dipole_length);
    buffer[i].swim.swimming = part.swimming().swimming;
#endif
    i++;
  }
}

void cuda_mpi_get_particles(
    const ParticleRange &particles,
    pinned_vector<CUDA_particle_data> &particle_data_host) {
  auto const n_part = particles.size();

  if (this_node > 0) {
    static std::vector<CUDA_particle_data> buffer;
    buffer.resize(n_part);
    /* pack local parts into buffer */
    pack_particles(particles, buffer.data());

    Utils::Mpi::gather_buffer(buffer, comm_cart);
  } else {
    particle_data_host.resize(n_part);

    /* Pack own particles */
    pack_particles(particles, particle_data_host.data());

    Utils::Mpi::gather_buffer(particle_data_host, comm_cart);
  }
}

/**
 * @brief Add a flat force (and torque) array to a range of particles.
 *
 * @param particles The particles the forces (and torques should be added to)
 * @param forces The forces as flat array of size 3 * particles.size()
 * @param torques The torques as flat array of size 3 * particles.size(),
 *                this is only touched if ROTATION is active.
 */
static void add_forces_and_torques(ParticleRange particles,
                                   Utils::Span<const float> forces,
                                   Utils::Span<const float> torques) {
  unsigned int i = 0;
  for (auto &p : particles) {
    for (unsigned int j = 0; j < 3; j++) {
      p.force()[j] += static_cast<double>(forces[3 * i + j]);
#ifdef ROTATION
      p.torque()[j] += static_cast<double>(torques[3 * i + j]);
#endif
    }
    i++;
  }
}

void cuda_mpi_send_forces(const ParticleRange &particles,
                          Utils::Span<float> host_forces,
                          Utils::Span<float> host_torques) {

  auto const size = 3 * particles.size();
  auto const n_elements = static_cast<int>(size);

  if (this_node > 0) {
    static std::vector<float> buffer_forces;
    static std::vector<float> buffer_torques;

    buffer_forces.resize(size);
    Utils::Mpi::scatter_buffer(buffer_forces.data(), n_elements, comm_cart);
#ifdef ROTATION
    buffer_torques.resize(size);
    Utils::Mpi::scatter_buffer(buffer_torques.data(), n_elements, comm_cart);
#endif
    add_forces_and_torques(particles, buffer_forces, buffer_torques);
  } else {
    Utils::Mpi::scatter_buffer(host_forces.data(), n_elements, comm_cart);
#ifdef ROTATION
    Utils::Mpi::scatter_buffer(host_torques.data(), n_elements, comm_cart);
#endif
    add_forces_and_torques(particles, host_forces, host_torques);
  }
}

void cuda_bcast_global_part_params() {
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
            sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
  EspressoSystemInterface::Instance().requestParticleStructGpu();
}

#endif // CUDA
