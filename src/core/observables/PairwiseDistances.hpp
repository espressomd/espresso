/*
 * Copyright (C) 2025 The ESPResSo project
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

#include "BoxGeometry.hpp"
#include "PidPairwiseDistancesObservable.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "particle_node.hpp"
#include "system/System.hpp"

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Observables {

/** @brief Track pairwise distances between two sets of particles. */
class PairwiseDistances : public PidPairwiseDistancesObservable {
public:
  using PidPairwiseDistancesObservable::PidPairwiseDistancesObservable;
  explicit PairwiseDistances(std::vector<int> const &ids,
                             std::vector<int> const &target_ids)
      : PidPairwiseDistancesObservable(ids, target_ids),
        m_pairs{get_unique_pairs(ids, target_ids)} {}

  /** @brief Evaluate the current contact times */
  std::vector<double>
  evaluate(boost::mpi::communicator const &comm, ParticleReferenceRange const &,
           ParticleObservables::traits<Particle> const &) const override {

    if (comm.rank() != 0) {
      return {};
    }
    // Get instances of the system and cell structures
    auto const &system = System::get_system();
    auto const &box_geo = *system.box_geo;
    auto &cell_structure = *system.cell_structure;

    std::vector<double> pairwise_distances;
    for (auto const &[pid1, pid2] : m_pairs) {
      auto const *p1 = cell_structure.get_local_particle(pid1);
      auto const *p2 = cell_structure.get_local_particle(pid2);
      auto const dist = box_geo.get_mi_vector(p1->pos(), p2->pos()).norm();
      pairwise_distances.emplace_back(dist);
    }
    return pairwise_distances;
  }

  std::vector<std::size_t> shape() const override { return {m_pairs.size()}; }

private:
  std::vector<std::pair<int, int>> m_pairs;

  std::vector<std::pair<int, int>>
  get_unique_pairs(std::vector<int> const &ids1, std::vector<int> const &ids2) {
    std::set<std::pair<int, int>> unique_pairs;
    for (int id1 : ids1) {
      for (int id2 : ids2) {
        if (id1 != id2) {
          unique_pairs.emplace(std::minmax(id1, id2));
        }
      }
    }
    return {unique_pairs.begin(), unique_pairs.end()};
  }
};

} // namespace Observables
