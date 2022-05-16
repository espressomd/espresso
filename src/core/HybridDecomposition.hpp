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

#ifndef ESPRESSO_HYBRID_DECOMPOSITION_HPP
#define ESPRESSO_HYBRID_DECOMPOSITION_HPP

#include "AtomDecomposition.hpp"
#include "Cell.hpp"
#include "Particle.hpp"
#include "ParticleDecomposition.hpp"
#include "RegularDecomposition.hpp"
#include "ghosts.hpp"
#include "integrate.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/optional.hpp>

#include <set>
#include <utility>
#include <vector>

/**
 * @brief Hybrid decomposition cell system.
 *
 * Store particles with short-range interactions
 * in a @ref RegularDecomposition cell system and
 * particles with long-range interactions
 * in a @ref AtomDecomposition (N-square) cell system.
 * All regular cells are coupled to the N-square cells.
 */
class HybridDecomposition : public ParticleDecomposition {
  boost::mpi::communicator m_comm;
  BoxGeometry m_box;
  double m_cutoff_regular;
  std::vector<Cell *> m_local_cells;
  std::vector<Cell *> m_ghost_cells;

  GhostCommunicator m_exchange_ghosts_comm;
  GhostCommunicator m_collect_ghost_force_comm;

  /** RegularDecomposition to hold the small particles */
  RegularDecomposition m_regular_decomposition;
  /** N-Square Decomposition to hold large particles */
  AtomDecomposition m_n_square;
  /** Set containing the types that should be handled using n_square */
  std::set<int> const m_n_square_types;

  bool is_n_square_type(int type_id) {
    return (m_n_square_types.find(type_id) != m_n_square_types.end());
  }

public:
  HybridDecomposition(boost::mpi::communicator comm, double cutoff_regular,
                      const BoxGeometry &box_geo,
                      const LocalBox<double> &local_box,
                      std::set<int> n_square_types);

  Utils::Vector3i get_cell_grid() const {
    return m_regular_decomposition.cell_grid;
  }

  Utils::Vector3d get_cell_size() const {
    return m_regular_decomposition.cell_size;
  }

  std::set<int> get_n_square_types() const { return m_n_square_types; }

  void resort(bool global, std::vector<ParticleChange> &diff) override;

  double get_cutoff_regular() const { return m_cutoff_regular; }

  GhostCommunicator const &exchange_ghosts_comm() const override {
    return m_exchange_ghosts_comm;
  }

  GhostCommunicator const &collect_ghost_force_comm() const override {
    return m_collect_ghost_force_comm;
  }

  Utils::Span<Cell *> local_cells() override {
    return Utils::make_span(m_local_cells);
  }

  Utils::Span<Cell *> ghost_cells() override {
    return Utils::make_span(m_ghost_cells);
  }

  Cell *particle_to_cell(Particle const &p) override {
    if (is_n_square_type(p.p.type)) {
      return m_n_square.particle_to_cell(p);
    }
    return m_regular_decomposition.particle_to_cell(p);
  }

  Utils::Vector3d max_cutoff() const override {
    return m_n_square.max_cutoff();
  }
  Utils::Vector3d max_range() const override { return m_n_square.max_range(); }

  boost::optional<BoxGeometry> minimum_image_distance() const override {
    return m_box;
  }

  BoxGeometry box() const override { return m_box; }

  Utils::Vector<std::size_t, 2> parts_per_decomposition_local() const;
};

#endif // ESPRESSO_HYBRID_DECOMPOSITION_HPP
