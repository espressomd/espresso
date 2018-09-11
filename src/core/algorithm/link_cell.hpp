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
#ifndef ALGORITHM_LINK_CELL_HPP
#define ALGORITHM_LINK_CELL_HPP

#include "utils/print.hpp"

namespace Algorithm {

/**
 * @brief Iterates over all particles in the cell range,
 *        and over all pairs within the cells and with
 *        their neighbors.
 */
template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename DistanceFunction>
void link_cell(CellIterator first, CellIterator last,
               ParticleKernel &&particle_kernel, PairKernel &&pair_kernel,
               DistanceFunction &&distance_function) {
  for (; first != last; ++first) {
    Utils::print("left", first->local_position[0], first->local_position[1],
                 first->local_position[2]);

    for (int i = 0; i != first->n; i++) {
      auto &p1 = first->part[i];

      particle_kernel(p1);

      /* Pairs in this cell */
      for (int j = i + 1; j < first->n; j++) {
        auto dist = distance_function(p1, first->part[j]);
        pair_kernel(p1, first->part[j], dist);
      }

      /* Pairs with neighbors */
      for (auto &neighbor : first->neighbors().red()) {
        Utils::print("right", neighbor->local_position[0],
                     neighbor->local_position[1], neighbor->local_position[2]);

        for (int j = 0; j < neighbor->n; j++) {
          auto &p2 = neighbor->part[j];
          auto dist = distance_function(p1, p2);
          pair_kernel(p1, p2, dist);
        }
      }
    }
  }
}
} // namespace Algorithm

#endif
