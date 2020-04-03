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

#include "cuda_interface.hpp"

#ifdef CUDA

#include "communication.hpp"
#include "config.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "serialization/CUDA_particle_data.hpp"

#include <utils/mpi/gather_buffer.hpp>
#include <utils/mpi/scatter_buffer.hpp>

/* TODO: We should only transfer data for enabled methods,
         not for those that are barely compiled in. (fw)
*/

static void pack_particles(ParticleRange particles,
                           CUDA_particle_data *buffer) {
  using Utils::Vector3f;

  int i = 0;
  for (auto const &part : particles) {
    buffer[i].p = static_cast<Vector3f>(folded_position(part.r.p, box_geo));

#ifdef CUDA
    buffer[i].identity = part.p.identity;
    buffer[i].v = static_cast<Vector3f>(part.m.v);
#ifdef VIRTUAL_SITES
    buffer[i].is_virtual = part.p.is_virtual;
#endif
#endif

#ifdef DIPOLES
    buffer[i].dip = static_cast<Vector3f>(part.calc_dip());
#endif

#if defined(LB_ELECTROHYDRODYNAMICS) && defined(CUDA)
    buffer[i].mu_E = static_cast<Vector3f>(part.p.mu_E);
#endif

#ifdef ELECTROSTATICS
    buffer[i].q = static_cast<float>(part.p.q);
#endif

#ifdef MASS
    buffer[i].mass = static_cast<float>(part.p.mass);
#endif

#ifdef ROTATION
    buffer[i].director = static_cast<Vector3f>(part.r.calc_director());
#endif

#ifdef ENGINE
    buffer[i].swim.v_swim = static_cast<float>(part.p.swim.v_swim);
    buffer[i].swim.f_swim = static_cast<float>(part.p.swim.f_swim);
    buffer[i].swim.director = buffer[i].director;

    buffer[i].swim.push_pull = part.p.swim.push_pull;
    buffer[i].swim.dipole_length =
        static_cast<float>(part.p.swim.dipole_length);
    buffer[i].swim.swimming = part.p.swim.swimming;
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

    Utils::Mpi::gather_buffer(buffer.data(), buffer.size(), comm_cart);
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
  int i = 0;
  for (auto &part : particles) {
    for (int j = 0; j < 3; j++) {
      part.f.f[j] += forces[3 * i + j];
#ifdef ROTATION
      part.f.torque[j] += torques[3 * i + j];
#endif
    }
    i++;
  }
}

void cuda_mpi_send_forces(const ParticleRange &particles,
                          Utils::Span<float> host_forces,
                          Utils::Span<float> host_torques) {
  auto const n_elements = 3 * particles.size();

  if (this_node > 0) {
    static std::vector<float> buffer_forces;
    static std::vector<float> buffer_torques;
    /* Alloc buffer */
    buffer_forces.resize(n_elements);

    Utils::Mpi::scatter_buffer(buffer_forces.data(), n_elements, comm_cart);
#ifdef ROTATION
    /* Alloc buffer */
    buffer_torques.resize(n_elements);

    Utils::Mpi::scatter_buffer(buffer_torques.data(), n_elements, comm_cart);
#endif

    add_forces_and_torques(particles, buffer_forces, buffer_torques);
  } else {
    /* Scatter forces to slaves */
    Utils::Mpi::scatter_buffer(host_forces.data(), n_elements, comm_cart);
#ifdef ROTATION
    Utils::Mpi::scatter_buffer(host_torques.data(), n_elements, comm_cart);
#endif

    add_forces_and_torques(particles, host_forces, host_torques);
  }
}

/** Takes a CUDA_energy struct and adds it to the core energy struct.
This cannot be done from inside cuda_common_cuda.cu:copy_energy_from_GPU()
because energy.hpp indirectly includes on mpi.h while .cu files may not depend
on mpi.h. */
void copy_CUDA_energy_to_energy(CUDA_energy energy_host) {
  energy.bonded[0] += energy_host.bonded;
  energy.non_bonded[0] += energy_host.non_bonded;
  if (energy.n_coulomb >= 1)
    energy.coulomb[0] += energy_host.coulomb;
  if (energy.n_dipolar >= 2)
    energy.dipolar[1] += energy_host.dipolar;
}

void cuda_bcast_global_part_params() { mpi_bcast_cuda_global_part_vars(); }

#endif /* ifdef CUDA */
