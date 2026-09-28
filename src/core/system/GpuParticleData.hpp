/*
 * Copyright (C) 2014-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_CUDA

#include "ParticleRange.hpp"
#include "ResourceCleanup.hpp"
#include "cuda/CudaHostAllocator.hpp"
#include "system/Leaf.hpp"

#include <utils/Vector.hpp>

#include <bitset>
#include <cstddef>
#include <memory>
#include <span>

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
 *
 * @note The AoS host-staging design (@ref pack_particles gathering fields
 * through the particle accessor API) remains the multi-rank CUDA data path.
 * The single-rank fast path fills per-field SoA staging buffers directly from
 * the @ref ParticleStore columns (@c pack_particles_soa), bypassing the AoS
 * pack + split kernel.
 */
class GpuParticleData : public System::Leaf<GpuParticleData>,
                        public std::enable_shared_from_this<GpuParticleData> {
public:
  /** @brief Particle properties that need to be communicated to the GPU. */
  struct prop {
    static constexpr std::size_t pos = 0;
    static constexpr std::size_t force = 1;
    static constexpr std::size_t torque = 2;
    static constexpr std::size_t q = 3;
    static constexpr std::size_t dip = 4;
    static constexpr std::size_t dip_fld = 5;
    using bitset = std::bitset<6>;
  };

  /** @brief Energies that are retrieved from the GPU. */
  struct GpuEnergy {
    float coulomb, dipolar;
  };

  /** @brief Subset of @ref Particle which is copied to the GPU. */
  struct GpuParticle {
    Utils::Vector3f p;
#ifdef ESPRESSO_DIPOLES
    Utils::Vector3f dip;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    float q;
#endif
    int identity;
  };

private:
  // forward declare
  class Storage;
  std::unique_ptr<Storage> m_data;
  void deinitialize() noexcept;
  using DeviceMemory =
      ResourceCleanup::Attorney<&GpuParticleData::deinitialize>;
  friend DeviceMemory;

  /** @brief Whether a device was found and data structures were allocated. */
  bool m_communication_enabled = false;
  /** @brief Whether to convert particle properties from AoS to SoA. */
  bool m_split_particle_struct = false;
  /** @brief Whether particle transfer to the GPU was requested. */
  bool m_need_particles_update = false;
  /** @brief Host and device containers. */

  bool has_compatible_device_impl() const;
  void gpu_init_particle_comm();
  void enable_particle_transfer();
  void copy_particles_to_device();
  void copy_particles_to_device(ParticleRange const &particles, int this_node,
                                bool single_rank);
  /** @brief Collect particles from all nodes to the head node. */
  void gather_particle_data(ParticleRange const &particles,
                            pinned_vector<GpuParticle> &particle_data_host,
                            int this_node);
  /**
   * @brief Fill the per-field SoA host staging buffers directly from the
   * @ref ParticleStore columns (single-rank fast path).
   *
   * Writes the enabled properties (position, charge, dipole moment) into the
   * @p positions / @p charges / @p dipoles spans in SoA layout, using the
   * same particle accessors and f64->f32 casts as @ref pack_particles feeding
   * @ref Storage::split_particle_struct, so the resulting device SoA buffers
   * are bit-identical to the AoS pack + split path. An empty span means the
   * property is not enabled and is skipped. Only meaningful on a single rank
   * (@c comm.size()==1): the store is rank-local, so the whole system is the
   * local range, and no MPI gather is needed.
   */
  void pack_particles_soa(ParticleRange const &particles,
                          std::span<float> positions, std::span<float> charges,
                          std::span<float> dipoles) const;
  void particles_scatter_forces(ParticleRange const &particles,
                                std::span<float> host_forces,
                                std::span<float> host_torques) const;
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  void particles_scatter_dip_fld(ParticleRange const &particles,
                                 std::span<float> host_dip_fld) const;
#endif

public:
  GpuParticleData();
  ~GpuParticleData();

  void update() {
    if (m_need_particles_update and m_communication_enabled) {
      copy_particles_to_device();
    }
  }
  void initialize();
  void enable_property(std::size_t property);
  void clear_energy_on_device();
  void copy_forces_to_host(ParticleRange const &particles, int this_node);
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  void copy_dip_fld_to_host(ParticleRange const &particles, int this_node);
#endif
  std::size_t n_particles() const;
  bool has_compatible_device() const;

  GpuEnergy copy_energy_to_host() const;
  GpuEnergy *get_energy_device() const;
  float *get_particle_positions_device() const;
  float *get_particle_forces_device() const;
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  float *get_particle_dip_fld_device() const;
#endif
#ifdef ESPRESSO_ROTATION
  float *get_particle_torques_device() const;
#endif
#ifdef ESPRESSO_DIPOLES
  float *get_particle_dipoles_device() const;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  float *get_particle_charges_device() const;
#endif
};

#endif // ESPRESSO_CUDA
