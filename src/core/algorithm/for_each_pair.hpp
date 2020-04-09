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
#ifndef CORE_ALGORITHM_PAIR_LOOP_HPP
#define CORE_ALGORITHM_PAIR_LOOP_HPP

#include <utility>

#include "link_cell.hpp"
#include "verlet_ia.hpp"

namespace Algorithm {
/**
 * @brief Run pair kernel for each particle pair from a cell range.
 *
 * Iterates over all cells in [first, last), and for every particle pair within
 * the cell and for each pair with the cells neighbors, @p distance_function is
 * evaluated and @p verlet_criterion is evaluated with the calculated distance.
 * Iff true, the @p pair_kernel is called.
 *
 * For details see @ref verlet_ia and @ref link_cell.
 *
 * Requirements on the types:
 * The Cell type has to provide a function %neighbors() that returns
 * a cell range comprised of the topological neighbors of the cell,
 * excluding the cell itself. The cells have to provide a %m_verlet_list
 * container that can be used to store particle pairs. It can be empty and is
 * not touched if @p use_verlet_list is false.
 *
 * <tt>verlet_criterion(p1, p2, distance_function(p1, p2))</tt> has to be valid
 * and convertible to bool.
 *
 * @tparam CellIterator     cell structure iterator
 * @tparam PairKernel       has to provide an %operator() member that can be
 *                          called with two particle references and a distance
 * @tparam DistanceFunction calculates distances between particles
 * @tparam VerletCriterion  determines whether a particle pair is interacting,
 *                          e.g. by comparing the particle pair distance to
 *                          the interaction range
 */
template <typename CellIterator, typename PairKernel, typename DistanceFunction,
          typename VerletCriterion>
void for_each_pair(CellIterator first, CellIterator last,
                   PairKernel &&pair_kernel,
                   DistanceFunction &&distance_function,
                   VerletCriterion &&verlet_criterion, bool use_verlet_list,
                   bool rebuild) {
  if (use_verlet_list) {
    verlet_ia(first, last, std::forward<PairKernel>(pair_kernel),
              std::forward<DistanceFunction>(distance_function),
              std::forward<VerletCriterion>(verlet_criterion), rebuild);
  } else {
    link_cell(first, last, std::forward<PairKernel>(pair_kernel),
              std::forward<DistanceFunction>(distance_function));
  }
}
} // namespace Algorithm

#endif
