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

#include <utility>

#include <boost/iterator/indirect_iterator.hpp>

#include "algorithm/for_each_pair.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "utils/Batch.hpp"

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  Vector3d vec21;
  double dist2;
};

namespace detail {
struct MinimalImageDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    Distance ret;
    get_mi_vector(ret.vec21, p1.r.p, p2.r.p);
    ret.dist2 = sqrlen(ret.vec21);

    return ret;
  }
};

struct LayeredMinimalImageDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    Distance ret;
    get_mi_vector(ret.vec21, p1.r.p, p2.r.p);
    ret.vec21[2] = p1.r.p[2] - p2.r.p[2];
    ret.dist2 = sqrlen(ret.vec21);

    return ret;
  }
};

struct EuclidianDistance {
  Distance operator()(Particle const &p1, Particle const &p2) const {
    Distance ret;
    ret.dist2 = distance2vec(p1.r.p, p2.r.p, ret.vec21);

    return ret;
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
  using Utils::make_batch;

  auto first = boost::make_indirect_iterator(local_cells.begin());
  auto last = boost::make_indirect_iterator(local_cells.end());

  /* In this case we reset l.p_old on the particles */
  if (rebuild_verletlist) {
    detail::decide_distance(
        first, last,
        /* Create a new functor that first runs the position
           copy and then the actual kernel. */
        make_batch([](Particle &p) { p.l.p_old = p.r.p; },
                   std::forward<ParticleKernel>(particle_kernel)),
        std::forward<PairKernel>(pair_kernel),
        VerletCriterion{skin, max_cut, coulomb_cutoff, dipolar_cutoff,
                        collision_detection_cutoff()});

    /* Now everything is up-to-date */
    rebuild_verletlist = 0;
  } else {
    detail::decide_distance(
        first, last, std::forward<ParticleKernel>(particle_kernel),
        std::forward<PairKernel>(pair_kernel),
        VerletCriterion{skin, max_cut, coulomb_cutoff, dipolar_cutoff});
  }
}

#endif
