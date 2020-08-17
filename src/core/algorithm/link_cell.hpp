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
#ifndef ALGORITHM_LINK_CELL_HPP
#define ALGORITHM_LINK_CELL_HPP

namespace Algorithm {

/**
 * @brief Iterates over all particles in the cell range,
 *        and over all pairs within the cells and with
 *        their neighbors.
 */
template <typename CellIterator, typename PairKernel>
void link_cell(CellIterator first, CellIterator last,
               PairKernel &&pair_kernel) {
  for (; first != last; ++first) {
    for (auto it = first->particles().begin(); it != first->particles().end();
         ++it) {
      auto &p1 = *it;

      /* Pairs in this cell */
      for (auto jt = std::next(it); jt != first->particles().end(); ++jt) {
        pair_kernel(p1, *jt);
      }

      /* Pairs with neighbors */
      for (auto &neighbor : first->neighbors().red()) {
        for (auto &p2 : neighbor->particles()) {
          pair_kernel(p1, p2);
        }
      }
    }
  }
}
} // namespace Algorithm

#endif
