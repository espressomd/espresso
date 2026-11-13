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

  /** @brief Number of bonds. */
  auto size() const { return m_bonds.size() - m_free_slots.size(); }
  bool empty() const { return size() == 0u; }

  void clear() {
    m_bonds.clear();
    m_free_slots.clear();
    m_slots_of_particle.clear();
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
        m_slots_of_particle[particle_ids[i]].push_back(slot);
      }
    }
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
    auto const it = m_slots_of_particle.find(particle_ids[0]);
    if (it == m_slots_of_particle.end()) {
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
    for (auto const slot : it->second) {
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
    auto const it = m_slots_of_particle.find(pid);
    if (it == m_slots_of_particle.end()) {
      return;
    }
    for (auto const slot : it->second) {
      kernel(std::as_const(m_bonds[static_cast<std::size_t>(slot)]));
    }
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
      auto const it = m_slots_of_particle.find(pid);
      if (it == m_slots_of_particle.end()) {
        continue; // particle listed several times in the bond
      }
      // keep the relative order of the remaining bonds of the particle
      std::erase(it->second, slot);
      if (it->second.empty()) {
        m_slots_of_particle.erase(it);
      }
    }
    bond.bond_id = -1;
    bond.particle_ids.clear();
    m_free_slots.push_back(slot);
  }

  /** @brief Bonds; slots are recycled and may be unused. */
  std::vector<Bond> m_bonds;
  /** @brief Unused slots of @ref m_bonds. */
  std::vector<int> m_free_slots;
  /** @brief For every particle, the slots of the bonds it participates in. */
  std::unordered_map<int, std::vector<int>> m_slots_of_particle;
};
