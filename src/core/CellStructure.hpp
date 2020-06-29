/*
 * Copyright (C) 2010-2020 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef ESPRESSO_CELLSTRUCTURE_HPP
#define ESPRESSO_CELLSTRUCTURE_HPP

#include "AtomDecomposition.hpp"
#include "BoxGeometry.hpp"
#include "Cell.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "ParticleDecomposition.hpp"
#include "ParticleList.hpp"
#include "ParticleRange.hpp"
#include "bond_error.hpp"
#include "ghosts.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/transform.hpp>

#include <vector>

/** Cell Structure */
enum CellStructureType : int {
  /** cell structure domain decomposition */
  CELL_STRUCTURE_DOMDEC = 1,
  /** cell structure n square */
  CELL_STRUCTURE_NSQUARE = 2
};

namespace Cells {
enum Resort : unsigned {
  RESORT_NONE = 0u,
  RESORT_LOCAL = 1u,
  RESORT_GLOBAL = 2u
};

/**
 * @brief Flags to select particle parts for communication.
 */
enum DataPart : unsigned {
  DATA_PART_NONE = 0u,       /**< Nothing */
  DATA_PART_PROPERTIES = 1u, /**< Particle::p */
  DATA_PART_POSITION = 2u,   /**< Particle::r */
  DATA_PART_MOMENTUM = 8u,   /**< Particle::m */
  DATA_PART_FORCE = 16u,     /**< Particle::f */
  DATA_PART_BONDS = 32u      /**< Particle::bonds */
};
} // namespace Cells

namespace Cells {
inline ParticleRange particles(Utils::Span<Cell *> cells) {
  /* Find first non-empty cell */
  auto first_non_empty =
      std::find_if(cells.begin(), cells.end(),
                   [](const Cell *c) { return not c->particles().empty(); });

  return {CellParticleIterator(first_non_empty, cells.end()),
          CellParticleIterator(cells.end())};
}
} // namespace Cells

/** Describes a cell structure / cell system. Contains information
 *  about the communication of cell contents (particles, ghosts, ...)
 *  between different nodes and the relation between particle
 *  positions and the cell system. All other properties of the cell
 *  system which are not common between different cell systems have to
 *  be stored in separate structures.
 */
struct CellStructure {
private:
  /** The local id-to-particle index */
  std::vector<Particle *> m_particle_index;
  /** Implementation of the primary particle decomposition */
  std::unique_ptr<ParticleDecomposition> m_decomposition =
      std::make_unique<AtomDecomposition>();
  /** Active type in m_decomposition */
  int m_type = CELL_STRUCTURE_NSQUARE;
  /** One of @ref Cells::Resort, announces the level of resort needed.
   */
  unsigned m_resort_particles = Cells::RESORT_NONE;

public:
  bool use_verlet_list = true;

  /**
   * @brief Update local particle index.
   *
   * Update the entry for a particle in the local particle
   * index.
   *
   * @param id Entry to update.
   * @param p Pointer to the particle.
   **/
  void update_particle_index(int id, Particle *p) {
    assert(id >= 0);
    assert(not p or id == p->identity());

    if (id >= m_particle_index.size())
      m_particle_index.resize(id + 1);

    m_particle_index[id] = p;
  }

  /**
   * @brief Update local particle index.
   *
   * Update the entry for a particle in the local particle
   * index.
   *
   * @param p Pointer to the particle.
   **/
  void update_particle_index(Particle &p) {
    update_particle_index(p.identity(), std::addressof(p));
  }

  /**
   * @brief Update local particle index.
   *
   * @param pl List of particles whose index entries should be updated.
   */
  void update_particle_index(ParticleList &pl) {
    for (auto &p : pl) {
      update_particle_index(p.identity(), std::addressof(p));
    }
  }

  /**
   * @brief Clear the particles index.
   */
  void clear_particle_index() { m_particle_index.clear(); }

private:
  /**
   * @brief Append a particle to a list and update this
   *        particle index accordingly.
   * @param pl List to add the particle to.
   * @param p Particle to add.
   */
  Particle &append_indexed_particle(ParticleList &pl, Particle &&p) {
    /* Check if cell may reallocate, in which case the index
     * entries for all particles in this cell have to be
     * updated. */
    auto const may_reallocate = pl.size() >= pl.capacity();
    auto &new_part = pl.insert(std::move(p));

    if (may_reallocate)
      update_particle_index(pl);
    else {
      update_particle_index(new_part);
    }

    return new_part;
  }

public:
  /**
   * @brief Get a local particle by id.
   *
   * @param id Particle to get.
   * @return Pointer to particle if it is local,
   *         nullptr otherwise.
   */
  Particle *get_local_particle(int id) {
    assert(id >= 0);

    if (id >= m_particle_index.size())
      return nullptr;

    return m_particle_index[id];
  }

  /** @overload */
  const Particle *get_local_particle(int id) const {
    assert(id >= 0);

    if (id >= m_particle_index.size())
      return nullptr;

    return m_particle_index[id];
  }

  template <class InputRange, class OutputIterator>
  void get_local_particles(InputRange ids, OutputIterator out) {
    boost::transform(ids, out,
                     [this](int id) { return get_local_particle(id); });
  }

public:
  int decomposition_type() const { return m_type; }

  /** Maximal pair range supported by current cell system. */
  Utils::Vector3d max_range() const;

  /** Return the global local_cells */
  Utils::Span<Cell *> local_cells();
  ParticleRange local_particles();
  ParticleRange ghost_particles();

  /** Cell system dependent function to find the right cell for a
   *  particle.
   *  \param  p Particle.
   *  \return pointer to cell where to put the particle, nullptr
   *          if the particle does not belong on this node.
   */
  Cell *particle_to_cell(const Particle &p);

  /**
   * @brief Add a particle.
   *
   * Moves a particle into the cell system. This adds
   * a particle to the local node, irrespective of where
   * it belongs.
   *
   * @param p Particle to add.
   * @return Pointer to the particle in the cell
   *         system.
   */
  Particle *add_particle(Particle &&p);

  /**
   * @brief Add a particle.
   *
   * Moves a particle into the cell system, if it
   * belongs to this node. Otherwise this does not
   * have an effect and the particle is discarded.
   * This can be used to add a particle without
   * knowledge where it should be placed by calling
   * the function on all nodes, it will then add
   * the particle in exactly one place.
   *
   * @param p Particle to add.
   * @return Pointer to particle if it is local, null
   *         otherwise.
   */
  Particle *add_local_particle(Particle &&p);

  /**
   * @brief Remove a particle.
   *
   * Removes a particle and all bonds pointing
   * to it. This is a collective call.
   *
   * @param id Id of particle to remove.
   */
  void remove_particle(int id);

  /**
   * @brief Get the maximal particle id on this node.
   *
   * This returns the highest particle id on
   * this node, or -1 if there are no particles on this node.
   */
  int get_max_local_particle_id() const;

  /**
   * @brief Remove all particles from the cell system.
   *
   * This allows linear time removal of all particles from
   * the system, removing each particle individually would
   * be quadratic.
   */
  void remove_all_particles();

  /**
   * @brief Get the underlyind particle decomposition.
   *
   * Should be used solely for informative purposes.
   *
   * @return The active particle decomposition.
   */
  const ParticleDecomposition &decomposition() const {
    return assert(m_decomposition), *m_decomposition;
  }

private:
  ParticleDecomposition &decomposition() {
    return assert(m_decomposition), *m_decomposition;
  }

public:
  /**
   * @brief Increase the local resort level at least to @p level.
   */
  void set_resort_particles(Cells::Resort level) {
    m_resort_particles |= level;
    assert(m_resort_particles >= level);
  }

  /**
   * @brief Get the currently scheduled resort level.
   */
  unsigned get_resort_particles() const { return m_resort_particles; }

  /**
   * @brief Set the resort level to sorted.
   */
  void clear_resort_particles() { m_resort_particles = Cells::RESORT_NONE; }

  /**
   * @brief Update ghost particles.
   *
   * This function updates the ghost particles with data
   * from the real particles.
   *
   * @param data_parts Particle parts to update, combination of @ref
   * Cells::DataPart
   */
  void ghosts_update(unsigned data_parts);
  /**
   * @brief Add forces from ghost particles to real particles.
   */
  void ghosts_reduce_forces();

private:
  /**
   * @brief Resolve ids to particles.
   *
   * @throws BondResolutionError if one of the ids
   *         was not found.
   *
   * @param partner_ids Ids to resolve.
   * @return Vector of Particle pointers.
   */
  auto resolve_bond_partners(Utils::Span<const int> partner_ids) {
    boost::container::static_vector<Particle *, 4> partners;
    get_local_particles(partner_ids, std::back_inserter(partners));

    /* Check if id resolution failed for any partner */
    if (boost::algorithm::any_of(
            partners, [](Particle *partner) { return partner == nullptr; })) {
      throw BondResolutionError{};
    }

    return partners;
  }

public:
  template <class Handler>
  void execute_bond_handler(Particle &p, Handler handler) {
    for (const BondView bond : p.bonds()) {
      auto const partner_ids = bond.partner_ids();

      try {
        auto partners = resolve_bond_partners(partner_ids);

        auto const bond_broken =
            handler(p, bond.bond_id(), Utils::make_span(partners));

        if (bond_broken) {
          bond_broken_error(p.identity(), partner_ids);
        }
      } catch (const BondResolutionError &) {
        bond_broken_error(p.identity(), partner_ids);
      }
    }
  }

private:
  /** Go through ghost cells and remove the ghost entries from the
      local particle index. */
  void invalidate_ghosts() {
    for (auto const &p : ghost_particles()) {
      if (get_local_particle(p.identity()) == &p) {
        update_particle_index(p.identity(), nullptr);
      }
    }
  }

public:
  /**
   * @brief Resort particles.
   */
  void resort_particles(int global_flag);

private:
  /** @brief Set the particle decomposition, keeping the particles. */
  void set_particle_decomposition(
      std::unique_ptr<ParticleDecomposition> &&decomposition) {
    clear_particle_index();

    auto local_parts = local_particles();
    std::vector<Particle> particles(local_parts.begin(), local_parts.end());

    m_decomposition = std::move(decomposition);

    for (auto &p : particles) {
      add_particle(std::move(p));
    }
  }

public:
  /**
   * @brief Set the particle decomposition to
   *        AtomDecomposition.
   *
   *        @param comm Communicator to use.
   *        @param box Box Geometry
   */
  void set_atom_decomposition(boost::mpi::communicator const &comm,
                              BoxGeometry const &box);

  /**
   * @brief Set the particle decomposition to
   *        DomainDecomposition.
   *
   *        @param comm Cartesian communicator to use.
   *        @param box Box Geometry
   *        @param local_geo Geoemtry of the local box.
   */
  void set_domain_decomposition(boost::mpi::communicator const &comm,
                                double range, BoxGeometry const &box,
                                LocalBox<double> const &local_geo);

  /**
   * @brief Return true if minimum image convention is
   *        needed for distance calculation. */
  bool minimum_image_distance() const {
    return m_decomposition->minimum_image_distance();
  }
};

#endif // ESPRESSO_CELLSTRUCTURE_HPP
