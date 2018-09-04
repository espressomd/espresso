/*
  Copyright (C) 2017-2018 The ESPResSo project

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

#ifndef CORE_UTILS_FOR_EACH_PAIR_HPP
#define CORE_UTILS_FOR_EACH_PAIR_HPP

#include <iterator>

namespace Utils {

/**
 * @brief Execute op for each pair of elements in (first, last] once.
 *
 * Diagonal elements are excluded. Pair are traversed ordered, so that
 * for op(*it, *jt), it holds that distance(it - first) < distance(jt - first),
 * and distance(it_n - first) < distance(it_n+1 - first) for consecutive calls.
 */
template <typename ForwardIterator, typename BinaryOp>
void for_each_pair(ForwardIterator first, ForwardIterator last, BinaryOp op) {
  while (first != last) {
    for (auto it = std::next(first); it != last; ++it) {
      op(*first, *it);
    }

    ++first;
  }
}

/*
 * @brief Range overload for for_each_pair.
 */
template <typename ForwardRange, typename BinaryOp>
void for_each_pair(ForwardRange &&rng, BinaryOp &&op) {
  using std::begin;
  using std::end;
  for_each_pair(begin(rng), end(rng), std::forward<BinaryOp>(op));
}
} // namespace Utils

#endif
