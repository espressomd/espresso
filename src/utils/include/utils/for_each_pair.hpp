/*
 * Copyright (C) 2017-2022 The ESPResSo project
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

#ifndef CORE_UTILS_FOR_EACH_PAIR_HPP
#define CORE_UTILS_FOR_EACH_PAIR_HPP

#include <iterator>

namespace Utils {

/**
 * @brief Execute op for each pair of elements in [first, last) once.
 *
 * Diagonal elements are excluded. Pairs are traversed ordered, so that
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

/** @overload */
template <typename ForwardRange, typename BinaryOp>
void for_each_pair(ForwardRange &&rng, BinaryOp &&op) {
  using std::begin;
  using std::end;
  for_each_pair(begin(rng), end(rng), std::forward<BinaryOp>(op));
}

/**
 * @brief Execute op for each pair of elements between [first1, last1) and
 * [first2, last2).
 *
 * Diagonal elements are *not* excluded. Pairs are traversed ordered, so that
 * for op(*it, *jt), it holds that distance(it - first) < distance(jt - first),
 * and distance(it_n - first) < distance(it_n+1 - first) for consecutive calls.
 */
template <typename ForwardIterator, typename BinaryOp>
void for_each_cartesian_pair(ForwardIterator first1, ForwardIterator last1,
                             ForwardIterator first2, ForwardIterator last2,
                             BinaryOp op) {
  while (first1 != last1) {
    for (auto it = first2; it != last2; ++it) {
      op(*first1, *it);
    }

    ++first1;
  }
}

/** @overload */
template <typename ForwardRange, typename BinaryOp>
void for_each_cartesian_pair(ForwardRange &&rng1, ForwardRange &&rng2,
                             BinaryOp &&op) {
  using std::begin;
  using std::end;
  for_each_cartesian_pair(begin(rng1), end(rng1), begin(rng2), end(rng2),
                          std::forward<BinaryOp>(op));
}

/**
 * @brief Execute op for each pair of elements between [first1, last1) and
 * [first2, last2) if a condition is satisfied.
 *
 * Diagonal elements are *not* excluded. Pairs are traversed ordered, so that
 * for op(*it, *jt), it holds that distance(it - first) < distance(jt - first),
 * and distance(it_n - first) < distance(it_n+1 - first) for consecutive calls.
 */
template <typename ForwardIterator, typename BinaryOp, typename BinaryCmp>
void for_each_cartesian_pair_if(ForwardIterator first1, ForwardIterator last1,
                                ForwardIterator first2, ForwardIterator last2,
                                BinaryOp op, BinaryCmp cmp) {
  while (first1 != last1) {
    for (auto it = first2; it != last2; ++it) {
      if (cmp(*first1, *it)) {
        op(*first1, *it);
      }
    }

    ++first1;
  }
}

/** @overload */
template <typename ForwardRange, typename BinaryOp, typename BinaryCmp>
void for_each_cartesian_pair_if(ForwardRange &&rng1, ForwardRange &&rng2,
                                BinaryOp &&op, BinaryCmp cmp) {
  using std::begin;
  using std::end;
  for_each_cartesian_pair_if(begin(rng1), end(rng1), begin(rng2), end(rng2),
                             std::forward<BinaryOp>(op),
                             std::forward<BinaryCmp>(cmp));
}
} // namespace Utils
#endif
