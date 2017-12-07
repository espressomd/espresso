/*
  Copyright (C) 2014,2015,2016 The ESPResSo project

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


#include "cuda_interface.hpp"

#ifdef CUDA

#include "utils/serialization/CUDA_particle_data.hpp"
#include "utils/mpi/gather_buffer.hpp"
#include "utils/mpi/scatter_buffer.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "debug.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "interaction_data.hpp"
#include "random.hpp"


/// MPI tag for cuda particle gathering
#define REQ_CUDAGETPARTS 0xcc01
/// MPI tag for cuda force gathering
#define REQ_CUDAGETFORCES 0xcc02

#ifdef ENGINE
static void cuda_mpi_send_v_cs_slave(ParticleRange particles);
#endif

void cuda_bcast_global_part_params() {
  COMM_TRACE(fprintf(stderr, "%d: cuda_bcast_global_part_params\n", this_node));
  mpi_bcast_cuda_global_part_vars();
  COMM_TRACE(fprintf(stderr, "%d: cuda_bcast_global_part_params finished\n",
                     this_node));
}

/* TODO: We should only transfer data for enabled methods,
         not for those that are barely compiled in. (fw)
*/

static void pack_particles(ParticleRange particles,
                           CUDA_particle_data *buffer) {
  int dummy[3] = {0, 0, 0};
  double pos[3];

  int i = 0;
  for (auto const &part : particles) {
    memmove(pos, part.r.p, 3 * sizeof(double));
    fold_position(pos, dummy);

    buffer[i].p[0] = static_cast<float>(pos[0]);
    buffer[i].p[1] = static_cast<float>(pos[1]);
    buffer[i].p[2] = static_cast<float>(pos[2]);

    buffer[i].v[0] = static_cast<float>(part.m.v[0]);
    buffer[i].v[1] = static_cast<float>(part.m.v[1]);
    buffer[i].v[2] = static_cast<float>(part.m.v[2]);
#ifdef IMMERSED_BOUNDARY
    buffer[i].isVirtual = part.p.isVirtual;
#endif

#ifdef DIPOLES
    buffer[i].dip[0] = static_cast<float>(part.r.dip[0]);
    buffer[i].dip[1] = static_cast<float>(part.r.dip[1]);
    buffer[i].dip[2] = static_cast<float>(part.r.dip[2]);
#endif

#ifdef SHANCHEN
    // SAW TODO: does this really need to be copied every time?
    for (int ii = 0; ii < 2 * LB_COMPONENTS; ii++) {
      buffer[i].solvation[ii] = static_cast<float>(part.p.solvation[ii]);
    }
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
    buffer[i].mu_E[0] = static_cast<float>(part.p.mu_E[0]);
    buffer[i].mu_E[1] = static_cast<float>(part.p.mu_E[1]);
    buffer[i].mu_E[2] = static_cast<float>(part.p.mu_E[2]);
#endif

#ifdef ELECTROSTATICS
    buffer[i].q = static_cast<float>(part.p.q);
#endif

#ifdef MASS
    buffer[i].mass = static_cast<float>(part.p.mass);
#endif

#ifdef ROTATION
    buffer[i].quatu[0] = static_cast<float>(part.r.quatu[0]);
    buffer[i].quatu[1] = static_cast<float>(part.r.quatu[1]);
    buffer[i].quatu[2] = static_cast<float>(part.r.quatu[2]);
#endif

#ifdef ENGINE
    buffer[i].swim.v_swim = static_cast<float>(part.swim.v_swim);
    buffer[i].swim.f_swim = static_cast<float>(part.swim.f_swim);
    buffer[i].swim.quatu[0] = static_cast<float>(part.r.quatu[0]);
    buffer[i].swim.quatu[1] = static_cast<float>(part.r.quatu[1]);
    buffer[i].swim.quatu[2] = static_cast<float>(part.r.quatu[2]);
#if defined(LB) || defined(LB_GPU)
    buffer[i].swim.push_pull = part.swim.push_pull;
    buffer[i].swim.dipole_length = static_cast<float>(part.swim.dipole_length);
#endif
    buffer[i].swim.swimming = part.swim.swimming;
#endif
    i++;
  }
}

void cuda_mpi_get_particles(ParticleRange particles,
                            CUDA_particle_data *particle_data_host) {
  auto const n_part = particles.size();

  if (this_node > 0) {
    COMM_TRACE(fprintf(stderr, "%d: get_particles_slave, %d particles\n",
                       this_node, n_part));
    static std::vector<CUDA_particle_data> buffer;
    buffer.resize(n_part);
    /* pack local parts into buffer */
    pack_particles(particles, buffer.data());

    Utils::Mpi::gather_buffer(buffer.data(), buffer.size(), comm_cart);
  } else {
    /* Pack own particles */
    pack_particles(particles, particle_data_host);

    Utils::Mpi::gather_buffer(particle_data_host, n_part, comm_cart);
  }

  COMM_TRACE(fprintf(stderr, "%d: finished get\n", this_node));
}

/**
 * @brief Add a flat force (and torque) array to a range of particles.
 *
 * @param particles The particles the forces (and torques should be added to)
 * @param forces The forces as flat array of size 3 * particls.size()
 * @param torques The torques as flat array of size 3 * particls.size(),
 *                this is only touched if ROTATION is active.
 */
static void add_forces_and_torques(ParticleRange particles, float *forces,
                                   float *torques) {
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

void cuda_mpi_send_forces(ParticleRange particles, float *host_forces,
                          float *host_torques = nullptr) {
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

    add_forces_and_torques(particles, buffer_forces.data(),
                           buffer_torques.data());
  } else {
    /* Scatter forces to slaves */
    Utils::Mpi::scatter_buffer(host_forces, n_elements, comm_cart);
#ifdef ROTATION
    Utils::Mpi::scatter_buffer(host_torques, n_elements, comm_cart);
#endif

    add_forces_and_torques(particles, host_forces, host_torques);
  }

  COMM_TRACE(fprintf(stderr, "%d: finished get\n", this_node));
}

#ifdef SHANCHEN
static void set_composition(ParticleRange particles,
                            CUDA_fluid_composition *composition) {
  for (auto &part : particles) {
    for (int ii = 0; ii < LB_COMPONENTS; ii++) {
      part.r.composition[ii] = composition->weight[ii];
    }
    composition++;
  }
}

void cuda_mpi_send_composition(ParticleRange particles,
                               CUDA_fluid_composition *host_composition) {
  auto const n_elements = particles.size();

  if (this_node > 0) {
    std::vector<CUDA_fluid_composition> buffer(n_elements);

    Utils::Mpi::scatter_buffer(buffer.data(), n_elements, comm_cart);
    set_composition(particles, buffer.data());
  } else {
    Utils::Mpi::scatter_buffer(host_composition, n_elements, comm_cart);
    set_composition(particles, host_composition);
  }
}
#endif /* SHANCHEN */

#if defined(ENGINE) && defined(LB_GPU)
namespace {
void set_v_cs(ParticleRange particles, CUDA_v_cs *v_cs) {
  for (auto &p : particles) {
    for (int i = 0; i < 3; i++) {
      p.swim.v_center[i] = v_cs->v_cs[0 + i];
      p.swim.v_source[i] = v_cs->v_cs[3 + i];
    }
    v_cs++;
  }
}
}

void cuda_mpi_send_v_cs(ParticleRange particles, CUDA_v_cs *host_v_cs) {
  // first collect number of particles on each node
  auto const n_part = particles.size();

  // call slave functions to provide the slave data
  if (this_node > 0) {
    std::vector<CUDA_v_cs> buffer(n_part);

    Utils::Mpi::scatter_buffer(buffer.data(), n_part, comm_cart);
    set_v_cs(particles, buffer.data());
  } else {
    Utils::Mpi::scatter_buffer(host_v_cs, n_part, comm_cart);
    set_v_cs(particles, host_v_cs);
  }
  COMM_TRACE(fprintf(stderr, "%d: finished send\n", this_node));
}
#endif // ifdef ENGINE

/** Takes a CUDA_energy struct and adds it to the core energy struct.
This cannot be done from inside cuda_common_cuda.cu:copy_energy_from_GPU()
because energy.hpp indirectly includes on mpi.h while .cu files may not depend
on mpi.h. */
void copy_CUDA_energy_to_energy(CUDA_energy energy_host) {
  energy.bonded[0] += energy_host.bonded;
  energy.non_bonded[0] += energy_host.non_bonded;
  if (energy.n_coulomb>=1)
    energy.coulomb[0] += energy_host.coulomb;
  if (energy.n_dipolar >=2)
    energy.dipolar[1] += energy_host.dipolar;
}

#endif /* ifdef CUDA */
