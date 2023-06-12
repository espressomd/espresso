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

#pragma once

#include "config/config.hpp"

#ifdef CUDA

#include "ParticleRange.hpp"
#include "cuda/CudaHostAllocator.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <bitset>
#include <cstddef>
#include <memory>

/**
 * @brief Particle data communication manager for the GPU.
 *
 * When data is synchronized between host and device memory, a subset
 * of the @ref Particle struct is copied from each particle on the host
 * to the corresponding @ref GpuParticle struct on the device via
 * @ref GpuParticleData::update(). Once the transfer is complete,
 * the particle AoS on the device is copied (or "split") to a SoA
 * automatically.
 *
 * Note that once a particle member is requested, the corresponding device
 * memory is allocated and populated at every time step, even when the GPU
 * method that originally requested the data is disabled.
 */
class GpuParticleData {
public:
  /** @brief Particle properties that need to be communicated to the GPU. */
  struct prop {
    static constexpr std::size_t pos = 0;
    static constexpr std::size_t force = 1;
    static constexpr std::size_t torque = 2;
    static constexpr std::size_t q = 3;
    static constexpr std::size_t dip = 4;
    using bitset = std::bitset<5>;
  };

  /** @brief Energies that are retrieved from the GPU. */
  struct GpuEnergy {
    float coulomb, dipolar;
  };

  /** @brief Subset of @ref Particle which is copied to the GPU. */
  struct GpuParticle {
    Utils::Vector3f p;
#ifdef DIPOLES
    Utils::Vector3f dip;
#endif
#ifdef ELECTROSTATICS
    float q;
#endif
    int identity;
  };

private:
  // forward declare
  class Storage;
  /** @brief Whether a device was found and data structures were allocated. */
  bool m_communication_enabled = false;
  /** @brief Whether to convert particle properties from AoS to SoA. */
  bool m_split_particle_struct = false;
  /** @brief Whether particle transfer to the GPU was requested. */
  bool m_need_particles_update = false;
  /** @brief Host and device containers. */
  std::shared_ptr<Storage> m_data;

  bool has_compatible_device_impl() const;
  void gpu_init_particle_comm();
  void enable_particle_transfer();
  void copy_particles_to_device();
  void copy_particles_to_device(ParticleRange const &particles, int this_node);
  /** @brief Collect particles from all nodes to the head node. */
  void gather_particle_data(ParticleRange const &particles,
                            pinned_vector<GpuParticle> &particle_data_host,
                            int this_node);
  void particles_scatter_forces(ParticleRange const &particles,
                                Utils::Span<float> host_forces,
                                Utils::Span<float> host_torques) const;

public:
  GpuParticleData() = default;
  ~GpuParticleData();

  void update() {
    if (m_need_particles_update and m_communication_enabled) {
      copy_particles_to_device();
    }
  }
  void init();
  void enable_property(std::size_t property);
  void clear_energy_on_device();
  void copy_forces_to_host(ParticleRange const &particles, int this_node);
  std::size_t n_particles() const;
  bool has_compatible_device() const;

  GpuEnergy copy_energy_to_host() const;
  GpuEnergy *get_energy_device() const;
  float *get_particle_positions_device() const;
  float *get_particle_forces_device() const;
#ifdef ROTATION
  float *get_particle_torques_device() const;
#endif
#ifdef DIPOLES
  float *get_particle_dipoles_device() const;
#endif
#ifdef ELECTROSTATICS
  float *get_particle_charges_device() const;
#endif
};

#endif // CUDA
