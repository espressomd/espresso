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
#ifndef CORE_ALGORITHM_VERLET_IA_HPP
#define CORE_ALGORITHM_VERLET_IA_HPP

#include <utility>

namespace Algorithm {
namespace detail {

template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename DistanceFunction, typename VerletCriterion>
void update_and_kernel(CellIterator first, CellIterator last,
                       ParticleKernel &&particle_kernel,
                       PairKernel &&pair_kernel,
                       DistanceFunction &&distance_function,
                       VerletCriterion &&verlet_criterion) {
  for (; first != last; ++first) {
    /* Clear the VL */
    first->m_verlet_list.clear();

    for (auto it = first->particles().begin(); it != first->particles().end();
         ++it) {
      auto &p1 = *it;

      particle_kernel(p1);

      /* Pairs in this cell */
      for (auto jt = std::next(it); jt != first->particles().end(); ++jt) {
        auto const dist = distance_function(p1, *jt);
        if (verlet_criterion(p1, *jt, dist)) {
          pair_kernel(p1, *jt, dist);
          first->m_verlet_list.emplace_back(&p1, &(*jt));
        }
      }

      /* Pairs with neighbors */
      for (auto &neighbor : first->neighbors().red()) {
        for (auto &p2 : neighbor->particles()) {
          auto dist = distance_function(p1, p2);
          if (verlet_criterion(p1, p2, dist)) {
            pair_kernel(p1, p2, dist);
            first->m_verlet_list.emplace_back(&p1, &p2);
          }
        }
      }
    }
  }
}

template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename DistanceFunction>
void kernel(CellIterator first, CellIterator last,
            ParticleKernel &&particle_kernel, PairKernel &&pair_kernel,
            DistanceFunction &&distance_function) {
  for (; first != last; ++first) {
    for (auto &p : first->particles()) {
      particle_kernel(p);
    }

    for (auto &pair : first->m_verlet_list) {
      auto const dist = distance_function(*pair.first, *pair.second);
      pair_kernel(*pair.first, *pair.second, dist);
    }
  }
}
} // namespace detail

/**
 * @brief Iterates over all particles in the cell range
 *        and all pairs in the Verlet list of the cells.
 *        If rebuild is true, all neighbor cells are iterated
 *        and the Verlet lists are updated with the so found pairs.
 */
template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename DistanceFunction, typename VerletCriterion>
void verlet_ia(CellIterator first, CellIterator last,
               ParticleKernel &&particle_kernel, PairKernel &&pair_kernel,
               DistanceFunction &&distance_function,
               VerletCriterion &&verlet_criterion, bool rebuild) {
  if (rebuild) {
    detail::update_and_kernel(first, last,
                              std::forward<ParticleKernel>(particle_kernel),
                              std::forward<PairKernel>(pair_kernel),
                              std::forward<DistanceFunction>(distance_function),
                              std::forward<VerletCriterion>(verlet_criterion));
  } else {
    detail::kernel(first, last, std::forward<ParticleKernel>(particle_kernel),
                   std::forward<PairKernel>(pair_kernel),
                   std::forward<DistanceFunction>(distance_function));
  }
}
} // namespace Algorithm

#endif
