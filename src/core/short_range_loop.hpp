/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef CORE_SHORT_RANGE_HPP
#define CORE_SHORT_RANGE_HPP

#include "algorithm/for_each_pair.hpp"
#include "cells.hpp"
#include "grid.hpp"

#include <boost/iterator/indirect_iterator.hpp>
#include <profiler/profiler.hpp>

#include <utility>

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  explicit Distance(Utils::Vector3d const &vec21)
      : vec21(vec21), dist2(vec21.norm2()) {}

  Utils::Vector3d vec21;
  double dist2;
};

namespace detail {
struct MinimalImageDistance {
  const BoxGeometry box;

  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(get_mi_vector(p1.r.p, p2.r.p, box));
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
        std::forward<PairKernel>(pair_kernel), MinimalImageDistance{box_geo},
        std::forward<VerletCriterion>(verlet_criterion),
        cell_structure.use_verlet_list, rebuild_verletlist);
    break;
  }
}

/**
 * @brief Functor that returns true for
 *        any arguments.
 */
struct True {
  template <class... T> bool operator()(T...) const { return true; }
};
} // namespace detail

template <class ParticleKernel, class PairKernel,
          class VerletCriterion = detail::True>
void short_range_loop(ParticleKernel &&particle_kernel,
                      PairKernel &&pair_kernel,
                      const VerletCriterion &verlet_criterion = {}) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (cell_structure.min_range != INACTIVE_CUTOFF) {
    auto first =
        boost::make_indirect_iterator(cell_structure.local_cells().begin());
    auto last =
        boost::make_indirect_iterator(cell_structure.local_cells().end());

    detail::decide_distance(
        first, last, std::forward<ParticleKernel>(particle_kernel),
        std::forward<PairKernel>(pair_kernel), verlet_criterion);

    rebuild_verletlist = false;
  } else {
    for (auto &p : cell_structure.local_particles()) {
      particle_kernel(p);
    }
  }
}

#endif
