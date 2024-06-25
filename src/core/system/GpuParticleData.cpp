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

#include "config/config.hpp"

#ifdef CUDA

#include "GpuParticleData.hpp"

#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "cuda/CudaHostAllocator.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>
#include <utils/mpi/scatter_buffer.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/split_free.hpp>

#include <cstddef>
#include <span>
#include <vector>

void GpuParticleData::enable_particle_transfer() {
  if (m_need_particles_update and not m_communication_enabled) {
    if (::this_node == 0) {
      gpu_init_particle_comm();
    }
    m_communication_enabled = true;
  }
}

void GpuParticleData::copy_particles_to_device() {
  auto const &cell_structure = *System::get_system().cell_structure;
  copy_particles_to_device(cell_structure.local_particles(), ::this_node);
}

bool GpuParticleData::has_compatible_device() const {
  auto result = false;
  if (::this_node == 0) {
    result = has_compatible_device_impl();
  }
  boost::mpi::broadcast(::comm_cart, result, 0);
  return result;
}

BOOST_IS_BITWISE_SERIALIZABLE(GpuParticleData::GpuParticle)

namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar, GpuParticleData::GpuParticle &p, unsigned const) {
  ar >> make_array(reinterpret_cast<char *>(&p),
                   sizeof(GpuParticleData::GpuParticle));
}
template <typename Archive>
void save(Archive &ar, GpuParticleData::GpuParticle const &p, unsigned const) {
  ar << make_array(reinterpret_cast<char const *>(&p),
                   sizeof(GpuParticleData::GpuParticle));
}
} // namespace serialization
} // namespace boost

BOOST_SERIALIZATION_SPLIT_FREE(GpuParticleData::GpuParticle)

static void pack_particles(ParticleRange const &particles,
                           GpuParticleData::GpuParticle *buffer) {
  auto const &box = *System::get_system().box_geo;
  std::size_t i = 0u;
  for (auto const &p : particles) {
    buffer[i].p = static_cast<Utils::Vector3f>(box.folded_position(p.pos()));
#ifdef DIPOLES
    buffer[i].dip = static_cast<Utils::Vector3f>(p.calc_dip());
#endif
#ifdef ELECTROSTATICS
    buffer[i].q = static_cast<float>(p.q());
#endif
    buffer[i].identity = p.id();
    i++;
  }
}

void GpuParticleData::gather_particle_data(
    ParticleRange const &particles,
    pinned_vector<GpuParticle> &particle_data_host, int this_node) {
  auto const n_part = particles.size();

  if (this_node > 0) {
    static std::vector<GpuParticle> buffer;
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
 * @param particles The particles the forces (and torques) should be added to
 * @param forces The forces as flat array of size 3 * particles.size()
 * @param torques The torques as flat array of size 3 * particles.size(),
 *                this is only touched if ROTATION is active.
 */
static void add_forces_and_torques(ParticleRange const &particles,
                                   std::span<const float> forces,
                                   std::span<const float> torques) {
  std::size_t i = 0ul;
  for (auto &p : particles) {
    for (std::size_t j = 0ul; j < 3ul; j++) {
      p.force()[j] += static_cast<double>(forces[3ul * i + j]);
#ifdef ROTATION
      p.torque()[j] += static_cast<double>(torques[3ul * i + j]);
#endif
    }
    i++;
  }
}

/**
 * @brief Distribute forces to the worker nodes, and add them to the particles.
 *
 * @param particles    The particles for which the forces (and torques) should
 *                     be added to.
 * @param host_forces  The forces as flat array of size 3 * particles.size(),
 *                     only relevant on the head node.
 * @param host_torques The torques as flat array of size 3 * particles.size(),
 *                     this is only touched if ROTATION is active. Only
 *                     relevant on the head node.
 */
void GpuParticleData::particles_scatter_forces(
    ParticleRange const &particles, std::span<float> host_forces,
    std::span<float> host_torques) const {

  auto const size = 3ul * particles.size();
  auto const n_elements = static_cast<int>(size);

  if (::this_node > 0) {
    static std::vector<float> buffer_forces;
    static std::vector<float> buffer_torques;

    buffer_forces.resize(size);
    Utils::Mpi::scatter_buffer(buffer_forces.data(), n_elements, ::comm_cart);
#ifdef ROTATION
    buffer_torques.resize(size);
    Utils::Mpi::scatter_buffer(buffer_torques.data(), n_elements, ::comm_cart);
#endif
    add_forces_and_torques(particles, buffer_forces, buffer_torques);
  } else {
    Utils::Mpi::scatter_buffer(host_forces.data(), n_elements, ::comm_cart);
#ifdef ROTATION
    Utils::Mpi::scatter_buffer(host_torques.data(), n_elements, ::comm_cart);
#endif
    add_forces_and_torques(particles, host_forces, host_torques);
  }
}

#endif
