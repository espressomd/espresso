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

#include "algorithm/link_cell.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "integrate.hpp"

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
 * @brief Functor that returns true for
 *        any arguments.
 */
struct True {
  template <class... T> bool operator()(T...) const { return true; }
};

template <class PairKernel, class DistanceFunction, class VerletCriterion>
void pair_loop(PairKernel &&pair_kernel, DistanceFunction df,
               const VerletCriterion &verlet_criterion) {
  auto first =
      boost::make_indirect_iterator(cell_structure.local_cells().begin());
  auto last = boost::make_indirect_iterator(cell_structure.local_cells().end());

  if (cell_structure.use_verlet_list && cell_structure.m_rebuild_verlet_list) {
    cell_structure.m_verlet_list.clear();

    Algorithm::link_cell(
        first, last,
        [&pair_kernel, &df, &verlet_criterion](Particle &p1, Particle &p2) {
          auto const d = df(p1, p2);
          if (verlet_criterion(p1, p2, d)) {
            cell_structure.m_verlet_list.emplace_back(&p1, &p2);
            pair_kernel(p1, p2, d);
          }
        });
    cell_structure.m_rebuild_verlet_list = false;
  } else if (cell_structure.use_verlet_list &&
             not cell_structure.m_rebuild_verlet_list) {
    for (auto &pair : cell_structure.m_verlet_list) {
      pair_kernel(*pair.first, *pair.second, df(*pair.first, *pair.second));
    }
  } else {
    Algorithm::link_cell(
        first, last,
        [&pair_kernel, &df, &verlet_criterion](Particle &p1, Particle &p2) {
          auto const d = df(p1, p2);
          if (verlet_criterion(p1, p2, d)) {
            cell_structure.m_verlet_list.emplace_back(&p1, &p2);
            pair_kernel(p1, p2, d);
          }
        });
  }
}

} // namespace detail

template <class ParticleKernel, class PairKernel,
          class VerletCriterion = detail::True>
void short_range_loop(ParticleKernel &&particle_kernel,
                      PairKernel &&pair_kernel,
                      const VerletCriterion &verlet_criterion = {}) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  auto local_particles = cell_structure.local_particles();
  for (auto &p : local_particles) {
    particle_kernel(p);
  }

  if (interaction_range() != INACTIVE_CUTOFF) {
    if (cell_structure.decomposition().minimum_image_distance()) {
      detail::pair_loop(std::forward<PairKernel>(pair_kernel),
                        detail::MinimalImageDistance{box_geo},
                        verlet_criterion);
    } else {
      detail::pair_loop(std::forward<PairKernel>(pair_kernel),
                        detail::EuclidianDistance{}, verlet_criterion);
    }
  }
}

#endif
