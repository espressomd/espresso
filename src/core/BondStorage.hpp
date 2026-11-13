/*
 * Copyright (C) 2026 The ESPResSo project
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

#include <boost/container/static_vector.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * @brief Global bond topology.
 *
 * Holds every bond of the system exactly once, independently of where the
 * participating particles are stored. Every MPI rank holds an identical
 * copy: all mutations must be executed on all ranks, with the same
 * arguments and in the same order (see @ref add_bond() and
 * @ref remove_bond()).
 *
 * A bond is a bond id (the key of the bonded interaction) and the ids of
 * the participating particles in creation order. The first particle is the
 * owner of the bond: the bond is evaluated on the MPI rank where the owner
 * resides. The same bond may be stored more than once.
 */
class BondStorage {
public:
  /** @brief Maximal number of particles in a bond. */
  static constexpr std::size_t max_particles = 4u;

  struct Bond {
    /** @brief Key of the bonded interaction, negative for unused slots. */
    int bond_id = -1;
    /** @brief Participants in creation order; the first one is the owner. */
    boost::container::static_vector<int, max_particles> particle_ids;

    int owner_id() const { return particle_ids.front(); }
    auto partner_ids() const {
      return std::span<const int>(particle_ids.data() + 1,
                                  particle_ids.size() - 1u);
    }
    bool has_participant(int pid) const {
      return std::ranges::find(particle_ids, pid) != particle_ids.end();
    }
  };

  /**
   * @brief Read-only view of the pair, angle and dihedral bonds, sorted by
   * owner.
   *
   * For each number of partners, the bonds are stored in one contiguous
   * table in compressed-row format, so that the bonds of a particle are
   * found without any search. Bonds without partners are not included.
   */
  class OwnedBonds {
  public:
    /** @brief Number of bonds owned by a particle, by number of partners. */
    struct Counts {
      int pair = 0;
      int angle = 0;
      int dihedral = 0;
    };

    /** @brief Number of pair, angle and dihedral bonds owned by a particle. */
    Counts counts(int pid) const {
      return {m_tables[0].count(pid), m_tables[1].count(pid),
              m_tables[2].count(pid)};
    }

    /**
     * @brief Execute a kernel for the bonds owned by a particle.
     *
     * Over all particles, every bond is visited exactly once.
     *
     * @param pid     Particle id.
     * @param kernel  Callable invoked as
     *                <tt>kernel(int bond_id, std::span<const int>
     *                partner_ids)</tt>, with the partners in creation order.
     */
    template <typename Kernel>
    void for_each_bond_of(int pid, Kernel &&kernel) const {
      for (std::size_t n_partners = 1u; n_partners <= 3u; ++n_partners) {
        auto const &table = m_tables[n_partners - 1u];
        auto const width = n_partners + 1u;
        auto const [begin, end] = table.range(pid);
        for (auto row = begin; row < end; ++row) {
          auto const *data = table.rows.data() + row * width;
          kernel(data[0], std::span<const int>(data + 1, n_partners));
        }
      }
    }

  private:
    friend class BondStorage;

    struct Table {
      /** @brief Rows of particle @c pid are
       *  <tt>[offsets[pid], offsets[pid + 1])</tt>. */
      std::vector<std::size_t> offsets;
      /** @brief Per bond, the bond id followed by the partner ids. */
      std::vector<int> rows;

      std::pair<std::size_t, std::size_t> range(int pid) const {
        auto const i = static_cast<std::size_t>(pid);
        if (pid < 0 or i + 1u >= offsets.size()) {
          return {0u, 0u};
        }
        return {offsets[i], offsets[i + 1u]};
      }
      int count(int pid) const {
        auto const [begin, end] = range(pid);
        return static_cast<int>(end - begin);
      }
    };

    void rebuild(std::vector<Bond> const &bonds) {
      int max_owner_id = -1;
      for (auto const &bond : bonds) {
        if (bond.bond_id >= 0) {
          max_owner_id = std::max(max_owner_id, bond.owner_id());
        }
      }
      auto const n_owners = static_cast<std::size_t>(max_owner_id + 1);
      for (std::size_t n_partners = 1u; n_partners <= 3u; ++n_partners) {
        auto &table = m_tables[n_partners - 1u];
        auto const width = n_partners + 1u;
        auto const selected = [width](Bond const &bond) {
          return bond.bond_id >= 0 and bond.particle_ids.size() == width;
        };
        // counting sort by owner
        table.offsets.assign(n_owners + 1u, 0u);
        for (auto const &bond : bonds) {
          if (selected(bond)) {
            ++table.offsets[static_cast<std::size_t>(bond.owner_id()) + 1u];
          }
        }
        for (std::size_t i = 0u; i < n_owners; ++i) {
          table.offsets[i + 1u] += table.offsets[i];
        }
        table.rows.resize(table.offsets[n_owners] * width);
        auto next_row = table.offsets;
        for (auto const &bond : bonds) {
          if (selected(bond)) {
            auto const owner = static_cast<std::size_t>(bond.owner_id());
            auto *data = table.rows.data() + next_row[owner]++ * width;
            data[0] = bond.bond_id;
            std::ranges::copy(bond.partner_ids(), data + 1);
          }
        }
      }
    }

    /** @brief Tables for the bonds with 1, 2 and 3 partners. */
    std::array<Table, 3> m_tables;
  };

  /** @brief Number of bonds. */
  auto size() const { return m_bonds.size() - m_free_slots.size(); }
  bool empty() const { return size() == 0u; }

  void clear() {
    m_bonds.clear();
    m_free_slots.clear();
    m_particles.clear();
    m_owned_bonds_dirty = true;
  }

  /**
   * @brief Add a bond.
   * @param bond_id       Key of the bonded interaction.
   * @param particle_ids  Participants in creation order, owner first.
   */
  void add(int bond_id, std::span<const int> particle_ids) {
    assert(bond_id >= 0);
    assert(not particle_ids.empty() and particle_ids.size() <= max_particles);
    int slot;
    if (m_free_slots.empty()) {
      slot = static_cast<int>(m_bonds.size());
      m_bonds.emplace_back();
    } else {
      slot = m_free_slots.back();
      m_free_slots.pop_back();
    }
    auto &bond = m_bonds[static_cast<std::size_t>(slot)];
    bond.bond_id = bond_id;
    bond.particle_ids.assign(particle_ids.begin(), particle_ids.end());
    for (std::size_t i = 0u; i < particle_ids.size(); ++i) {
      // a particle listed several times in a bond is indexed once
      auto const first = std::ranges::find(particle_ids, particle_ids[i]);
      if (first == particle_ids.begin() + static_cast<std::ptrdiff_t>(i)) {
        m_particles[particle_ids[i]].slots.push_back(slot);
      }
    }
    m_owned_bonds_dirty = true;
  }

  /**
   * @brief Remove one bond.
   *
   * A bond matches when it has the given bond id and the same participants.
   * A bond that is also owned by @p particle_ids[0] is preferred, which
   * disambiguates bonds stored twice with swapped ownership; the order of
   * the other participants never matters.
   *
   * @return whether a bond was removed.
   */
  bool remove(int bond_id, std::span<const int> particle_ids) {
    if (particle_ids.empty() or particle_ids.size() > max_particles) {
      return false;
    }
    auto const it = m_particles.find(particle_ids[0]);
    if (it == m_particles.end()) {
      return false;
    }
    using Ids = boost::container::static_vector<int, max_particles>;
    auto const sorted = [](std::span<const int> ids) {
      Ids out(ids.begin(), ids.end());
      std::ranges::sort(out);
      return out;
    };
    auto const wanted = sorted(particle_ids);
    int fallback = -1;
    int match = -1;
    for (auto const slot : it->second.slots) {
      auto const &bond = m_bonds[static_cast<std::size_t>(slot)];
      if (bond.bond_id != bond_id or sorted(bond.particle_ids) != wanted) {
        continue;
      }
      if (bond.owner_id() == particle_ids[0]) {
        match = slot;
        break;
      }
      if (fallback == -1) {
        fallback = slot;
      }
    }
    if (match == -1) {
      match = fallback;
    }
    if (match == -1) {
      return false;
    }
    erase_slot(match);
    return true;
  }

  /**
   * @brief Execute a kernel for all bonds a particle participates in,
   * whether it owns them or not.
   *
   * The storage must not be modified by the kernel.
   *
   * @param pid     Particle id.
   * @param kernel  Callable invoked as <tt>kernel(Bond const &)</tt>.
   */
  template <typename Kernel>
  void for_each_bond_of(int pid, Kernel &&kernel) const {
    auto const it = m_particles.find(pid);
    if (it == m_particles.end()) {
      return;
    }
    for (auto const slot : it->second.slots) {
      kernel(std::as_const(m_bonds[static_cast<std::size_t>(slot)]));
    }
  }

  /**
   * @brief The bonds sorted by owner, for fast traversal.
   *
   * Rebuilt here if the storage was modified since the last call, with a
   * cost proportional to the number of bonds; this function must therefore
   * not be called concurrently. The returned object is safe to read
   * concurrently, and is valid until the storage is modified.
   */
  OwnedBonds const &owned_bonds() const {
    if (m_owned_bonds_dirty) {
      m_owned_bonds.rebuild(m_bonds);
      m_owned_bonds_dirty = false;
    }
    return m_owned_bonds;
  }

  /**
   * @brief Execute a kernel for all bonds.
   *
   * The storage must not be modified by the kernel.
   *
   * @param kernel  Callable invoked as <tt>kernel(Bond const &)</tt>.
   */
  template <typename Kernel> void for_each_bond(Kernel &&kernel) const {
    for (auto const &bond : m_bonds) {
      if (bond.bond_id >= 0) {
        kernel(bond);
      }
    }
  }

private:
  void erase_slot(int slot) {
    auto &bond = m_bonds[static_cast<std::size_t>(slot)];
    for (auto const pid : bond.particle_ids) {
      auto const it = m_particles.find(pid);
      if (it == m_particles.end()) {
        continue; // particle listed several times in the bond
      }
      // keep the relative order of the remaining bonds of the particle
      std::erase(it->second.slots, slot);
      if (it->second.slots.empty()) {
        m_particles.erase(it);
      }
    }
    bond.bond_id = -1;
    bond.particle_ids.clear();
    m_free_slots.push_back(slot);
    m_owned_bonds_dirty = true;
  }

  /** @brief Bonds; slots are recycled and may be unused. */
  std::vector<Bond> m_bonds;
  /** @brief Unused slots of @ref m_bonds. */
  std::vector<int> m_free_slots;
  struct ParticleEntry {
    /** @brief Slots of the bonds the particle participates in. */
    std::vector<int> slots;
  };

  /** @brief Bonds of every particle that participates in a bond. */
  std::unordered_map<int, ParticleEntry> m_particles;
  /** @brief Cache of the bonds sorted by owner, see @ref owned_bonds(). */
  mutable OwnedBonds m_owned_bonds;
  mutable bool m_owned_bonds_dirty = false;
};
