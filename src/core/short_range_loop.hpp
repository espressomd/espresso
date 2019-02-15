/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CORE_SHORT_RANGE_HPP
#define CORE_SHORT_RANGE_HPP

#include "algorithm/for_each_pair.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <boost/iterator/indirect_iterator.hpp>
#include <profiler/profiler.hpp>

#include <utility>

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  explicit Distance(Vector3d const &vec21)
      : vec21(vec21), dist2(vec21.norm2()) {}

  Vector3d vec21;
  const double dist2;
};

namespace detail {
struct MinimalImageDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(get_mi_vector(p1.r.p, p2.r.p));
  }
};

struct LayeredMinimalImageDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    auto mi_dist = get_mi_vector(p1.r.p, p2.r.p);
    mi_dist[2] = p1.r.p[2] - p2.r.p[2];

    return Distance(mi_dist);
  }
};

struct EuclidianDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(p1.r.p - p2.r.p);
  }
};

/**
 * @brief Decided which distance function to use depending on the
          cell system, and call the pair code.
*/
template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename VerletCriterion>
void decide_distance(CellIterator first, CellIterator last,
                     ParticleKernel &&particle_kernel, PairKernel &&pair_kernel,
                     VerletCriterion &&verlet_criterion) {
  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    Algorithm::for_each_pair(
        first, last, std::forward<ParticleKernel>(particle_kernel),
        std::forward<PairKernel>(pair_kernel), EuclidianDistance{},
        std::forward<VerletCriterion>(verlet_criterion),
        cell_structure.use_verlet_list, rebuild_verletlist);
    break;
  case CELL_STRUCTURE_NSQUARE:
    Algorithm::for_each_pair(
        first, last, std::forward<ParticleKernel>(particle_kernel),
        std::forward<PairKernel>(pair_kernel), MinimalImageDistance{},
        std::forward<VerletCriterion>(verlet_criterion),
        cell_structure.use_verlet_list, rebuild_verletlist);
    break;
  case CELL_STRUCTURE_LAYERED:
    Algorithm::for_each_pair(
        first, last, std::forward<ParticleKernel>(particle_kernel),
        std::forward<PairKernel>(pair_kernel), LayeredMinimalImageDistance{},
        std::forward<VerletCriterion>(verlet_criterion),
        cell_structure.use_verlet_list, rebuild_verletlist);
    break;
  }
}
} // namespace detail

template <typename ParticleKernel, typename PairKernel>
void short_range_loop(ParticleKernel &&particle_kernel,
                      PairKernel &&pair_kernel) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  auto first = boost::make_indirect_iterator(local_cells.begin());
  auto last = boost::make_indirect_iterator(local_cells.end());

  detail::decide_distance(
      first, last, std::forward<ParticleKernel>(particle_kernel),
      std::forward<PairKernel>(pair_kernel),
      VerletCriterion{skin, max_cut, coulomb_cutoff, dipolar_cutoff,
                      collision_detection_cutoff()});

  rebuild_verletlist = 0;
}

#endif
