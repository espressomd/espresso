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

#include "Cell.hpp"
#include "Particle.hpp"
#include "ParticleList.hpp"
#include "ParticleRange.hpp"
#include "ghosts.hpp"

#include <boost/range/algorithm/find_if.hpp>
#include <vector>

/** Cell Structure */
enum {
  /** Flag indicating that there is no cell system yet. Only at the
   *  VERY beginning of the code startup.
   */
  CELL_STRUCTURE_NONEYET = -1,
  /** Flag indicating that the current cell structure will be used further on */
  CELL_STRUCTURE_CURRENT = 0,
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
}

/** List of cell pointers. */
struct CellPList {
  ParticleRange particles() const {
    /* Find first non-empty cell */
    auto first = std::find_if(cell, cell + n,
                              [](const Cell *c) { return not c->empty(); });

    return {CellParticleIterator(first, cell + n, 0),
            CellParticleIterator(cell + n, cell + n, 0)};
  }

  Cell **begin() { return cell; }
  Cell **end() { return cell + n; }

  Cell *operator[](int i) { return assert(i < n), cell[i]; }

  Cell **cell = nullptr;
  int n = 0;
};

/** Describes a cell structure / cell system. Contains information
 *  about the communication of cell contents (particles, ghosts, ...)
 *  between different nodes and the relation between particle
 *  positions and the cell system. All other properties of the cell
 *  system which are not common between different cell systems have to
 *  be stored in separate structures.
 */
struct CellStructure {
private:
  std::vector<Particle *> m_particle_index;

public:
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
   * @brief Update local particle index.
   *
   * @param pl List of particles whose index entries should be updated.
   */
  void update_particle_index(ParticleList *pl) {
    assert(pl), update_particle_index(*pl);
  }

private:
  /**
   * @brief Append a particle to a list and update this
   *        particle index accordingly.
   * @param pl List to add the particle to.
   * @param p Particle to add.
   */
  void append_indexed_particle(ParticleList *pl, Particle &&p) {
    auto const old_data = pl->data();
    pl->push_back(std::move(p));

    /* If the list storage moved due to reallocation,
     * we have to update the index for all particles,
     * otherwise just for the particle that we added. */
    if (old_data != pl->data())
      update_particle_index(pl);
    else {
      update_particle_index(pl->back());
    }
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

  /** @overload get_local_particle */
  const Particle *get_local_particle(int id) const {
    assert(id >= 0);

    if (id >= m_particle_index.size())
      return nullptr;

    return m_particle_index[id];
  }

  std::vector<Cell *> m_local_cells = {};
  std::vector<Cell *> m_ghost_cells = {};

  /** type descriptor */
  int type = CELL_STRUCTURE_NONEYET;

  bool use_verlet_list = true;

  /** Maximal pair range supported by current cell system. */
  Utils::Vector3d max_range = {};

  /** Minimum range that has to be supported. */
  double min_range;

  /** Return the global local_cells */
  CellPList local_cells() {
    return {m_local_cells.data(), static_cast<int>(m_local_cells.size())};
  }
  /** Return the global ghost_cells */
  CellPList ghost_cells() {
    return {m_ghost_cells.data(), static_cast<int>(m_ghost_cells.size())};
  }

  /** Communicator to exchange ghost particles. */
  GhostCommunicator exchange_ghosts_comm;
  /** Communicator to collect ghost forces. */
  GhostCommunicator collect_ghost_force_comm;

  /** Cell system dependent function to find the right cell for a
   *  particle.
   *  \param  p Particle.
   *  \return pointer to cell where to put the particle, nullptr
   *          if the particle does not belong on this node.
   */
  Cell *(*particle_to_cell)(const Particle &p) = nullptr;

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
   * to it. This is a colective call.
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

private:
  /** One of @ref Cells::Resort, announces the level of resort needed.
   */
  unsigned m_resort_particles = Cells::RESORT_NONE;

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
};

#endif // ESPRESSO_CELLSTRUCTURE_HPP
