/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "particle_store/MigrationPack.hpp"
#include "particle_store/ParticleStore.hpp"
#include "system/System.hpp"

#include <boost/mpi/communicator.hpp>

#include <array>
#include <optional>
#include <utility>
#include <vector>

/**
 * @brief A particle VIEW that owns the standalone ParticleStore backing its
 * single row.
 *
 * A particle copied to the head node lands in a single-row standalone store
 * (populated via @ref MigrationPack::unpack_rows, the per-field wire the
 * production head-node fetch uses); this class bundles that store so its
 * lifetime is tied to the copy, which matters because the store holds Kokkos
 * Views that must be destroyed before Kokkos::finalize(). Mimics the subset of
 * the std::optional<Particle> interface used by the tests (bool, has_value,
 * operator*).
 */
class ParticleWithStore {
  ParticleStore m_store{};
  std::optional<Particle> m_particle{};

  void reattach() {
    if (m_particle) {
      // Re-point the view at the moved store's row 0; the values already live
      // in the moved-from store's columns/sidecars.
      m_particle->attach_to_store(m_store, 0);
    }
  }

public:
  ParticleWithStore() = default;
  ParticleWithStore(ParticleWithStore &&other) noexcept
      : m_store(std::move(other.m_store)),
        m_particle(std::move(other.m_particle)) {
    // the moved store keeps its Kokkos Views; re-point the view at it
    reattach();
  }
  ParticleWithStore &operator=(ParticleWithStore &&other) noexcept {
    m_store = std::move(other.m_store);
    m_particle = std::move(other.m_particle);
    reattach();
    return *this;
  }
  ParticleWithStore(ParticleWithStore const &) = delete;
  ParticleWithStore &operator=(ParticleWithStore const &) = delete;

  /** Populate the single-row store from a per-field wire buffer produced by
   *  @ref MigrationPack::pack_rows for a single row, and bind a view to it. */
  void unpack_from(std::vector<char> const &buffer) {
    m_store.begin_rebuild(1u, 0u);
    m_store.finish_rebuild();
    MigrationPack::unpack_rows(m_store, 0, buffer);
    m_particle = m_store.make_view(0);
  }

  /** Copy a live source row (same-rank fast path) into the single-row store and
   *  bind a view to it. */
  void copy_from(ParticleStore const &source, int const source_row) {
    m_store.begin_rebuild(1u, 0u);
    m_store.finish_rebuild();
    m_store.copy_row(source, source_row, 0);
    m_particle = m_store.make_view(0);
  }

  explicit operator bool() const { return m_particle.has_value(); }
  bool has_value() const { return m_particle.has_value(); }
  Particle &operator*() { return *m_particle; }
  Particle const &operator*() const { return *m_particle; }
  Particle *operator->() { return &*m_particle; }
  Particle const *operator->() const { return &*m_particle; }
};

/**
 * @brief Copy a particle (with every per-particle field) to the head node.
 *
 * All fields live in the @ref ParticleStore columns and are transmitted via
 * the per-field @ref MigrationPack wire (the same wire the production
 * head-node fetch cache uses). The returned copy owns a standalone single-row
 * store so its accessors work.
 */
inline ParticleWithStore
copy_particle_to_head_node(boost::mpi::communicator const &comm,
                           System::System &system, int p_id) {
  ParticleWithStore result{};
  system.cell_structure->ensure_particle_store_synchronized();
  auto p = system.cell_structure->get_local_particle(p_id);
  auto const owns_it = p and not p->is_ghost();

  if (comm.rank() == 0) {
    if (owns_it) {
      // Same-rank fast path: copy the live row directly.
      result.copy_from(*p->store(), p->store_row());
    } else {
      // Receive the per-field migration buffer from the owning rank.
      std::vector<char> buffer;
      comm.recv(boost::mpi::any_source, p_id, buffer);
      result.unpack_from(buffer);
    }
  } else if (owns_it) {
    // Owning worker rank: pack the single row and send it to the head node.
    std::array<int, 1> const rows{p->store_row()};
    std::vector<char> buffer;
    MigrationPack::pack_rows(*p->store(), rows, buffer);
    comm.send(0, p_id, buffer);
  }
  return result;
}
