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
 *
 * The bonds are stored per owner, in a table indexed by particle id, so
 * that the bonds a particle owns are found without any search. A bond is
 * found from its other participants through links to its owner.
 */
class BondStorage {
public:
  /** @brief Maximal number of particles in a bond. */
  static constexpr std::size_t max_particles = 4u;

  struct Bond {
    /** @brief Key of the bonded interaction. */
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
   * @brief Read-only view of the pair, angle and dihedral bonds of each
   * owner.
   *
   * Bonds without partners are not included. Always up to date with the
   * storage, and safe to read concurrently as long as the storage is not
   * modified.
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
      auto const *owner = find(pid);
      return (owner) ? owner->counts : Counts{};
    }

    /**
     * @brief Execute a kernel for the bonds owned by a particle, in
     * insertion order.
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
      if (auto const *owner = find(pid)) {
        for (auto const &bond : owner->bonds) {
          if (bond.particle_ids.size() > 1u) {
            kernel(bond.bond_id, bond.partner_ids());
          }
        }
      }
    }

  private:
    friend class BondStorage;

    struct Owner {
      /** @brief Bonds owned by the particle, in insertion order. */
      std::vector<Bond> bonds;
      /** @brief Number of bonds in @ref bonds, by number of partners. */
      Counts counts;

      void update_counts(Bond const &bond, int delta) {
        switch (bond.particle_ids.size()) {
        case 2u:
          counts.pair += delta;
          break;
        case 3u:
          counts.angle += delta;
          break;
        case 4u:
          counts.dihedral += delta;
          break;
        default:
          break;
        }
      }
    };

    Owner const *find(int pid) const {
      auto const i = static_cast<std::size_t>(pid);
      return (pid >= 0 and i < m_owners.size()) ? &m_owners[i] : nullptr;
    }
    Owner *find(int pid) {
      return const_cast<Owner *>(std::as_const(*this).find(pid));
    }

    /** @brief Bonds of every particle, indexed by particle id. */
    std::vector<Owner> m_owners;
  };

  /** @brief Number of bonds. */
  auto size() const { return m_size; }
  bool empty() const { return size() == 0u; }

  void clear() {
    m_owned_bonds.m_owners.clear();
    m_links.clear();
    m_size = 0u;
  }

  /**
   * @brief Add a bond.
   * @param bond_id       Key of the bonded interaction.
   * @param particle_ids  Participants in creation order, owner first.
   */
  void add(int bond_id, std::span<const int> particle_ids) {
    assert(bond_id >= 0);
    assert(not particle_ids.empty() and particle_ids.size() <= max_particles);
    assert(std::ranges::all_of(particle_ids, [](int pid) { return pid >= 0; }));
    auto &owners = m_owned_bonds.m_owners;
    auto const owner_id = particle_ids[0];
    if (static_cast<std::size_t>(owner_id) >= owners.size()) {
      owners.resize(static_cast<std::size_t>(owner_id) + 1u);
    }
    auto &owner = owners[static_cast<std::size_t>(owner_id)];
    auto &bond = owner.bonds.emplace_back();
    bond.bond_id = bond_id;
    bond.particle_ids.assign(particle_ids.begin(), particle_ids.end());
    owner.update_counts(bond, +1);
    for_each_linked_partner(bond, [&](int pid) {
      auto &links = m_links[pid];
      auto const it = std::ranges::find(links, owner_id, &Link::owner_id);
      if (it != links.end()) {
        ++it->n_bonds;
      } else {
        links.push_back({owner_id, 1});
      }
    });
    ++m_size;
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
    using Ids = boost::container::static_vector<int, max_particles>;
    auto const sorted = [](std::span<const int> ids) {
      Ids out(ids.begin(), ids.end());
      std::ranges::sort(out);
      return out;
    };
    auto const wanted = sorted(particle_ids);
    auto const remove_from = [&](OwnedBonds::Owner &owner) {
      auto const it = std::ranges::find_if(owner.bonds, [&](Bond const &bond) {
        return bond.bond_id == bond_id and
               bond.particle_ids.size() == wanted.size() and
               sorted(bond.particle_ids) == wanted;
      });
      if (it == owner.bonds.end()) {
        return false;
      }
      erase_bond(owner, it);
      return true;
    };
    auto const pid = particle_ids[0];
    if (auto *owner = m_owned_bonds.find(pid)) {
      if (remove_from(*owner)) {
        return true;
      }
    }
    if (auto const it = m_links.find(pid); it != m_links.end()) {
      for (auto const &link : it->second) {
        // the links are modified by a removal: leave the loop right after
        if (remove_from(*m_owned_bonds.find(link.owner_id))) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * @brief Execute a kernel for all bonds a particle participates in,
   * whether it owns them or not.
   *
   * The bonds the particle owns are visited first, in insertion order.
   * The storage must not be modified by the kernel.
   *
   * @param pid     Particle id.
   * @param kernel  Callable invoked as <tt>kernel(Bond const &)</tt>.
   */
  template <typename Kernel>
  void for_each_bond_of(int pid, Kernel &&kernel) const {
    any_bond_of(pid, [&](Bond const &bond) {
      kernel(bond);
      return false;
    });
  }

  /**
   * @brief Check whether a particle participates in a bond accepted by
   * @p pred. Safe to call concurrently.
   *
   * @param pid   Particle id.
   * @param pred  Callable <tt>bool(Bond const &)</tt>.
   */
  template <typename Pred> bool any_bond_of(int pid, Pred &&pred) const {
    if (auto const *owner = m_owned_bonds.find(pid)) {
      if (std::ranges::any_of(owner->bonds, pred)) {
        return true;
      }
    }
    auto const it = m_links.find(pid);
    if (it == m_links.end()) {
      return false;
    }
    return std::ranges::any_of(it->second, [&](Link const &link) {
      auto const &bonds = m_owned_bonds.find(link.owner_id)->bonds;
      return std::ranges::any_of(bonds, [&](Bond const &bond) {
        return bond.has_participant(pid) and pred(bond);
      });
    });
  }

  /**
   * @brief The bonds of each owner, for fast traversal.
   *
   * The returned object is safe to read concurrently, as long as the
   * storage is not modified.
   */
  OwnedBonds const &owned_bonds() const { return m_owned_bonds; }

  /**
   * @brief Execute a kernel for all bonds.
   *
   * The storage must not be modified by the kernel.
   *
   * @param kernel  Callable invoked as <tt>kernel(Bond const &)</tt>.
   */
  template <typename Kernel> void for_each_bond(Kernel &&kernel) const {
    for (auto const &owner : m_owned_bonds.m_owners) {
      for (auto const &bond : owner.bonds) {
        kernel(bond);
      }
    }
  }

private:
  /** @brief Bonds of an owner in which another particle participates. */
  struct Link {
    int owner_id;
    /** @brief Number of such bonds. */
    int n_bonds;
  };

  /**
   * @brief Execute a kernel for the participants of a bond that reach it
   * through a link, i.e. once for every particle other than the owner.
   */
  template <typename Kernel>
  static void for_each_linked_partner(Bond const &bond, Kernel &&kernel) {
    auto const partner_ids = bond.partner_ids();
    for (auto it = partner_ids.begin(); it != partner_ids.end(); ++it) {
      if (*it != bond.owner_id() and
          std::find(partner_ids.begin(), it, *it) == it) {
        kernel(*it);
      }
    }
  }

  void erase_bond(OwnedBonds::Owner &owner,
                  std::vector<Bond>::iterator bond_it) {
    auto const bond = *bond_it;
    // keep the relative order of the remaining bonds of the owner
    owner.bonds.erase(bond_it);
    owner.update_counts(bond, -1);
    for_each_linked_partner(bond, [&](int pid) {
      auto const map_it = m_links.find(pid);
      auto &links = map_it->second;
      auto const it =
          std::ranges::find(links, bond.owner_id(), &Link::owner_id);
      if (--it->n_bonds == 0) {
        links.erase(it);
        if (links.empty()) {
          m_links.erase(map_it);
        }
      }
    });
    --m_size;
  }

  /** @brief Bonds, stored with their owner. */
  OwnedBonds m_owned_bonds;
  /** @brief For every particle, the owners of the bonds it participates in
   *  without owning them. */
  std::unordered_map<int, std::vector<Link>> m_links;
  std::size_t m_size = 0u;
};
